function output = createRJMCMC(profilesOutput,priorsPropOutput,varargin)

%% default parameters
iterations = 5000 ;
lsquaresType = 'wls' ;

proposalsType = 'uniform' ;
proposalsNumSelRegsPropLambda = [] ; 
proposalsMaxDimJumpPr = 0.2 ;
proposalsWDimJumRatio = 5/6 ;
proposalsRwSd = 10 ;
proposalsDeg_rwSd = 0.05 ;
proposalsDeg_sampled = 1 ;
mcmcPlotsOn = 0 ;

%% process inputs
if mod(length(varargin),2)
    error('Syntax error on command line. Enter Name Value pairs');
end
for i=1:2:length(varargin)-1
    switch varargin{i}
        case ('proposalsType')
            proposalsType = varargin{i+1};
        case('proposalsNumSelRegsPropLambda')
            proposalsNumSelRegsPropLambda = varargin{i+1};
        case('proposalsMaxDimJumpPr')
            proposalsMaxDimJumpPr = varargin{i+1};
        case('proposalsWDimJumRatio')
            proposalsWDimJumRatio = varargin{i+1};
        case('proposalsRwSd')
            proposalsRwSd = varargin{i+1};
        case ('proposalsDeg_rwSd')
            proposalsDeg_rwSd = varargin{i+1};
        case ('proposalsDeg_sampled')
            proposalsDeg_sampled = varargin{i+1};  
        case ('iterations')
            iterations = varargin{i+1};
        case('lsquaresType')
            lsquaresType = varargin{i+1}; 
        case('mcmcPlotsOn') 
            mcmcPlotsOn = varargin{i+1}; 
        otherwise
            error(['Unknown property ' varargin{i}  ' in the command line']);
    end
end

%% get various inputs 
regulatorNames = profilesOutput.regulatorNames ;
experiments = profilesOutput.experiments ;

numRegsPrior = priorsPropOutput.priors.nregulatorsPriorOutput.probs ;
regPrior = { priorsPropOutput.priors.regsPriorOutput.type , priorsPropOutput.priors.regsPriorOutput.regPriorsMatrixAllcomb} ;
thresholdPrior = { priorsPropOutput.priors.thresholdsPriorOutput.type , priorsPropOutput.priors.thresholdsPriorOutput.threshPriorDensity2} ;

%%% numbers
numExps = size(experiments,1);
numRegs = length(regulatorNames);
% check if any reg is rejected from priors
regCheck=ones(numRegs,1);
if strcmp(thresholdPrior{1},'gradient')
    for z=1:numRegs
        if regPrior{2}{1}(z)==0
            regCheck(z)=0;
        end
        if sum(thresholdPrior{2}(:,z))==0
            regCheck(z)=0;
        end
    end
end
allRegs = find(regCheck) ;
numRegs=length(allRegs);

origWarnings = warning;
warning('off','all');


%%% size of functions created via smoothing etc.
timeFuncSize = profilesOutput.timeFuncSize;
exprFuncSize = priorsPropOutput.profileFuncSize;

%%% various hyperparameters
%%% deg
deg = priorsPropOutput.priors.deg ; % initial value for deg
deg_sampled = proposalsDeg_sampled ;
minDegRate = priorsPropOutput.priors.minDegRate ; % minimum value for the bounded proposal ration for deg
maxDegRate = priorsPropOutput.priors.maxDegRate ; % maximum value for the bounded proposal ration for deg
prior_deg_a0 = priorsPropOutput.priors.prior_deg_a0 ; %
prior_deg_b0 = priorsPropOutput.priors.prior_deg_b0 ; % prior hyperparameters for gamma prior of deg
deg_rwSd = proposalsDeg_rwSd ;

%%% precision
precDeg0 = priorsPropOutput.priors.precDeg0 ;
precSigma0 = priorsPropOutput.priors.precSigma0 ;


%% get proposals distr
% mean (lambda) of the poisson proposal for the move probs
maxNumRegs = priorsPropOutput.maxNumRegs ;
maxNumRegs = min(maxNumRegs,numRegs);

if isempty(proposalsNumSelRegsPropLambda)
    numSelRegsPropLambda = round(maxNumRegs/2) ;
end
%%% move probs
maxDimJumpPr = proposalsMaxDimJumpPr ;

r = 1:maxNumRegs-1 ;
c = maxDimJumpPr/2-0.05:0.001:maxDimJumpPr ;
lc = length(c) ;
b = zeros(lc,maxNumRegs-1) ;
d = zeros(lc,maxNumRegs-1) ;
bd = zeros(lc,maxNumRegs -1) ;

for i=1:lc
    ci = c(i) ;
    for j=1:maxNumRegs-1
        rj = r(j) ;
        b(i,j) = ci*min(1,truncPoisson(rj+1,numSelRegsPropLambda,maxNumRegs)/truncPoisson(rj,numSelRegsPropLambda,maxNumRegs)) ;
        d(i,j) = ci*min(1,truncPoisson(rj-1,numSelRegsPropLambda,maxNumRegs)/truncPoisson(rj,numSelRegsPropLambda,maxNumRegs)) ;
        bd(i,j) = b(i,j) + d(i,j) ;
    end
end

if (maxNumRegs-1 > 1)
    idx = find(sum(bd' >= maxDimJumpPr),1,'first')-1 ;
    c = c(idx) ;
    changeDimConst = c ;
else
    idx = find(bd <= maxDimJumpPr,1,'last') ;
    c = c(idx) ;
    changeDimConst = c ;
end

sr = zeros(maxNumRegs+1,1) ;
tr = zeros(maxNumRegs+1,1) ;
br = zeros(maxNumRegs+1,1) ;
dr = zeros(maxNumRegs+1,1) ;
moveProbs = zeros(maxNumRegs+1,4) ;

sr(1) = 0 ;
tr(1) = 0 ;
br(1) = 1 ;
dr(1) = 0 ;
moveProbs(1,:) = cumsum([tr(1) sr(1) br(1) dr(1)]) ;

p = proposalsWDimJumRatio ;
if (maxNumRegs == numRegs)
    dr(end) = maxDimJumpPr ;
    br(end) = 0 ;
    tr(end) = 1 - dr(end) ;
    sr(end) = 0 ;
    moveProbs(end,:) = cumsum([tr(end) sr(end) br(end) dr(end)]) ;
elseif (maxNumRegs < numRegs)
    dr(end) = maxDimJumpPr ;
    br(end) = 0 ;
    tr(end) = p*(1 - maxDimJumpPr) ;
    sr(end) = 1 - maxDimJumpPr - tr(end) ;
    moveProbs(end,:) = cumsum([tr(end) sr(end) br(end) dr(end)]) ;
elseif (maxNumRegs > numRegs)
    dr(numRegs+1:end) = maxDimJumpPr ;
    br(numRegs+1:end) = 0 ;
    tr(numRegs+1:end) = 1 - maxDimJumpPr ;
    sr(numRegs+1:end) = 0 ;
    for jj=numRegs+1:maxNumRegs+1
        moveProbs(jj,:) = cumsum([tr(jj) sr(jj) br(jj) dr(jj)]) ;
    end
end

for j=2:maxNumRegs
    r = j-1 ;
    br(j) = changeDimConst*min(1,truncPoisson(r+1,numSelRegsPropLambda,maxNumRegs)/truncPoisson(r,numSelRegsPropLambda,maxNumRegs)) ;
    dr(j) = changeDimConst*min(1,truncPoisson(r-1,numSelRegsPropLambda,maxNumRegs)/truncPoisson(r,numSelRegsPropLambda,maxNumRegs)) ;
    bdr = br(j) + dr(j) ;
    tr(j) = p*(1 - bdr) ;
    sr(j) = 1 - bdr - tr(j) ;
    moveProbs(j,:) = cumsum([tr(j) sr(j) br(j) dr(j)]) ;
end
pr = moveProbs ;


%%% collect experiments' & other variables
ntimepointsAll = profilesOutput.ntimepointsAll ;
replicatesAll = profilesOutput.nreplicatesAll ;
timepointsAll = profilesOutput.timepointsAll ;
timeAll = profilesOutput.timeAll ;
childDataAll = profilesOutput.targetDataAll ;
eHatMatAll = profilesOutput.eHatMatAll ;
sampleSizeAll = profilesOutput.sampleSizeAll ;

%%% precision parameters
precPar_a1 = zeros(numExps,1) ;
for x = 1:numExps
    precPar_a1(x) = (precDeg0 + sampleSizeAll(x))/2 ;
end
precPar_b1 = zeros(numExps,1) ;

e_hat = cell(numExps,1) ;
prec = zeros(numExps,1) ;
for x=1:numExps
    e_hat{x} = ones(replicatesAll(x)*ntimepointsAll(x),1) ;
    prec(x) = 1 ;
    
end
sse = zeros(numExps,1) ;

%%% regulators expression, activation, etc.
regProfileAll = profilesOutput.regProfileAll ;
regProfileActivationAll = profilesOutput.regProfileActivationAll ;
regProfileActivationAllProposed = profilesOutput.regProfileActivationAllProposed ;

%%% min time between switches
minSwTime = profilesOutput.minSwTime ;

%%% delay
hSteps =  profilesOutput.hSteps ; % the delay time

%%% min/max's of regProfile of each reg
minRegProfiles =  profilesOutput.minRegProfiles ;
maxRegProfiles = profilesOutput.maxRegProfiles  ;

%%% std of random walk for threshold
rwSd = proposalsRwSd ;

%%% Y matrix for child data
YAll = profilesOutput.YAll;

%%% get the initial value of the regProfile threshold
regProfiles = cell2mat(regProfileAll) ;
regProfileThreshold = nanmean(regProfiles) ;   % initial threshold
regProfileMean = regProfileThreshold ; %  mean of the normal prior

%%% least squares parameters
lsquares.type = lsquaresType ;
if (strcmp(lsquares.type,'wls'))
    lsquares.weightsPhi = zeros(numExps,1) ;
    lsquares.weights = cell(numExps,1) ;
    lsquares.Q = cell(numExps,1) ;
    for x=1:numExps
        experiment = experiments{x,1};
        timepoints = experiment.ntimepoints; 
        w = 1./profilesOutput.targetSmoothData{x} ; 
        w = w/(prod(w)^(1/timepoints)) ;
        lsquares.weights{x} = w ;
    end
    lsquares.QAll = nan ;
    lsquaresWeightsPhiHistory = nan(iterations,numExps) ;
end

%% create variables for log-likelihood and chains
lLikelihoods = zeros(numExps,1);
lLikelihoodsProposed = zeros(numExps,1);
lLikelihoodHistory = nan(iterations,numExps+1);
likelihoodRatioHistory = nan(iterations,1);
degHistory = nan(iterations,1);
regProfileThresholdHistory = nan(iterations,numRegs);
precHistory = nan(iterations,numExps);
m0History = nan(iterations,numExps);
switchesHistory = cell(iterations,numExps);
statesHistory = cell(iterations,numExps);
birthratesHistory = cell(iterations,1);

%% space for switches and states
switchesAll = cell(numExps,1);
statesAll = cell(numExps,1);

targetmRNA = cell(numExps,1);

%% regSet initial parameters
selRegs = [] ; %  allRegs ; %
numSelRegs = length(selRegs) ;
%%% prior for initially selRegs
if (strcmp(regPrior{1},'range') && numSelRegs>0 )
    regPriorProbs = regPrior{2} ;
    [~,idxCombiSelRegs,~] = intersect(nchoosek(allRegs,numSelRegs),selRegs,'rows') ;
    priorSelRegs = regPriorProbs{numSelRegs}(idxCombiSelRegs,1) ;
else
    priorSelRegs = 1 ;
end
nonSelRegs = setdiff(allRegs,selRegs) ;
numNonSelRegs = length(nonSelRegs) ;
selRegsHistory = cell(iterations,1) ;
numSelRegsHistory = zeros(iterations,1) ;
selParentHistoryPc = zeros(numRegs,1) ;
numSelRegsHistoryPc = zeros(numRegs,1) ;

%% number of accepted moves
% M, S, B, D, deg respectively
accept = zeros(5,1);

%%% number of times a move of each kind is proposed
% M, S, B, D respectively
proposalNum = zeros(4,1);
    
accRatio = 0 ;

moveIdxHistory = zeros(iterations,1) ;

%%% prior for thresholds of initially selRegs
if strcmp(thresholdPrior{1},'gradient')
    levels = zeros(exprFuncSize,length(regCheck)) ;
    priorThresholds = zeros(length(regCheck),1) ;
    for zz = 1:numRegs
        z=allRegs(zz);
        levels(:,z) = linspace(minRegProfiles(z),maxRegProfiles(z),exprFuncSize);
        idxTh = find(regProfileThreshold(z) >= levels(:,z),1,'last') ;
        priorThresholds(z) = thresholdPrior{2}(idxTh,z) ;
    end
end


%%%%%%%  sampler starts %%%%%%%%%%%%%%%%%%%%%
for i = 1:iterations
        
    %%% get current model likelihood
    for x = 1:numExps
        regProfileSelRegs = regProfileAll{x}(:,selRegs);
        regProfileActivationSelRegs = regProfileActivationAll{x}(:,selRegs);
        time = timeAll{x};
        
        for g = 1:numSelRegs
            
            % get selected reg
            z = selRegs(g) ;
            
            %%% get activation
            regProfileActivationSelRegs(hSteps(x):end,g) = getRegProfileActivation(regProfileSelRegs(hSteps(x):end,g),regProfileThreshold(z));
            % set the act func before delay as the first after delay to avoid switches
            regProfileActivationSelRegs(1:hSteps(x)-1,g) = regProfileActivationSelRegs(hSteps(x),g) ;
        end
        
        regProfileActivationAll{x}(:,selRegs) = regProfileActivationSelRegs;
        
        %%% get vector with state at each (dense) time point for all regs
        twosMat = repmat(pow2(0:(numSelRegs-1)),timeFuncSize,1);
        stateChanges = getGeneActivation(regProfileActivationSelRegs,twosMat);
        %%% get states used and switch times
        %%% if too fast switches, get them to the largest switch and set
        %%% the state as after the last
        [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x)); % states used and at which times
        switchesAll{x} = switches;
        statesAll{x} = states;
    end
    %%% get the states used across experiments
    statesUsed = getStatesUsed(statesAll);
    %     numStatesProposed = length(statesUsed) ;
    
    %%%% perfrom multi-exp regression
    % get m0's and tau-rates
    qAll = [] ;
    if (strcmp(lsquares.type,'wls'))        
        for x=1:numExps
            yphi = log(prec(x)*e_hat{x}.^2) ;
            xphi = -repmat(log(lsquares.weights{x}),replicatesAll(x),1) ;
            phi_hat = ((xphi'*xphi)\xphi')*yphi ;
            phi_hat = max(min(phi_hat,1),0);
            lsquares.weightsPhi(x) = phi_hat ;
            qx = repmat(lsquares.weights{x}.^(-lsquares.weightsPhi(x)),replicatesAll(x),1) ;
            qAll = [qAll ; qx]   ;
            lsquares.Q{x} = diag(qx) ;            
        end
        QAll = diag(qAll) ;
        lsquares.QAll = QAll ;        
    end
    [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,replicatesAll,lsquares) ;
    m0All(m0All<0)=0;
    birthRates(birthRates<0)=0;
    
    for x = 1:numExps
        
        br = zeros(length(statesAll{x}),1) ;
        for y=1:length(br)
            br(y) = find(statesAll{x}(y) == statesUsed) ;
        end
        targetmRNA{x} = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
        
        %%% calculate SSE of profile to data
        e_hat{x} = calculateEHat(targetmRNA{x},childDataAll{x},eHatMatAll{x},replicatesAll(x)) ;
        
        if (strcmp(lsquares.type,'wls'))
            QAllx = lsquares.Q{x} ;
            sse(x) = e_hat{x}'*(QAllx\e_hat{x}) ;
        else
            sse(x) = sum(e_hat{x} .^ 2) ;
        end
        
        precPar_b1(x) = 2/(sampleSizeAll(x)*sse(x)/(sampleSizeAll(x)-length(br)-1) + precDeg0*precSigma0^2) ;
        
        %calculate precision (used for sigma_sq) - Gibbs
        prec(x) = gamrnd(precPar_a1(x),precPar_b1(x)) ;
        
        if strcmp(lsquares.type,'wls')
            lLikelihoods(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
        else
            lLikelihoods(x) = logLikelihood(e_hat{x},prec(x))  ;
        end
        
    end
    lLikelihood = sum(lLikelihoods) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% choose one of the possible moves %%%%%%%%%%%%%%
    
    %%% first compute the probabilities of each move given the current
    %%% number of selected regs
    p_r = pr(numSelRegs+1,:) ;
    
    %%% choose which move
    u = rand ; % 4 ; %
    moveIdx = find(p_r >= u,1,'first') ;
    
    %%% do the M (move Threshold) jump
    if (moveIdx == 1)
        
        %%% add another M jump
        proposalNum(1) = proposalNum(1) + 1 ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% try changing regProfile threshold %%%%%%%%%%
        %%% pick one of the selected regs and change its threshold
        selParentIdx = randi(numSelRegs) ; %
        regIdx = selRegs(selParentIdx) ;
        
        regProfileThresholdProposed = regProfileThreshold;
        
        th = fastnormrnd(regProfileThreshold(regIdx),rwSd,1) ;
        while (th <= minRegProfiles(regIdx) || th >= maxRegProfiles(regIdx))
            th = fastnormrnd(regProfileThreshold(regIdx),rwSd,1) ;
        end
        regProfileThresholdProposed(regIdx) = th ;
        
        %%% this uses a prior vaguely centered in the middle of the
        %%% expression levels.
        priorRatioNumSelRegs = 1 ;
        priorRatioSelRegs = 1 ;
        if (strcmp(thresholdPrior{1},'gradient'))
            priorThresholdsProposed = priorThresholds ;
            idxTh = find(th >= levels(:,regIdx),1,'last' ) ;
            priorThresholdsProposed(regIdx) = thresholdPrior{2}(idxTh,regIdx) ;            
            priorRatioThresholds = priorThresholdsProposed(regIdx)/priorThresholds(regIdx) ;
        else
            priorRatioThresholds = 1 ;
        end
        priorRatio = priorRatioNumSelRegs*priorRatioSelRegs*priorRatioThresholds ;
        proposalRatio = 1 ;
        
        if (priorRatio ~= 0 && ~isnan(priorRatio) && ~isinf(priorRatio))
            

            for x = 1:numExps
                
                time = timeAll{x};
                
                % get the activation function of all the selected regs
                regProfileActivationSelRegsProposed = regProfileActivationAll{x}(:,selRegs);
                
                % change the activation function of the reg with moved threshold
                % first get the expression of the reg to have its threshold moved
                regProfileMoveParent = regProfileAll{x}(:,regIdx);
                
                % use the function to get its activation
                regProfileActivationSelRegsProposed(hSteps(x):end,selParentIdx) = getRegProfileActivation(regProfileMoveParent(hSteps(x):end),th);
                % set the act func before delay as the first after delay to avoid switches
                regProfileActivationSelRegsProposed(1:hSteps(x)-1,selParentIdx) = regProfileActivationSelRegsProposed(hSteps(x),selParentIdx) ;
                
                % change the activation function of this reg
                regProfileActivationAllProposed{x}(:,regIdx) = regProfileActivationSelRegsProposed(:,selParentIdx);
                
                % get the states of the selected regs
                twosMat = repmat(pow2(0:(numSelRegs-1)),timeFuncSize,1);
                stateChanges = getGeneActivation(regProfileActivationSelRegsProposed,twosMat);
                [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x)) ;
                switchesAll{x} = switches;
                statesAll{x} = states;
            end
            
            statesUsed = getStatesUsed(statesAll) ;
            
            [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,replicatesAll,lsquares);
            m0All(m0All<0)=0;
            birthRates(birthRates<0)=0;
            
            
            for x = 1:numExps
                br = zeros(length(statesAll{x}),1) ;
                for y=1:length(br)
                    br(y) = find(statesAll{x}(y) == statesUsed) ;
                end
                targetmRNAProposed = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                
                e_hat{x} = calculateEHat(targetmRNAProposed,childDataAll{x},eHatMatAll{x},replicatesAll(x));
                
                if strcmp(lsquares.type,'wls')
                    lLikelihoodsProposed(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
                else
                    lLikelihoodsProposed(x) = logLikelihood(e_hat{x},prec(x))  ;
                end
            end
            lLikelihoodProposed = sum(lLikelihoodsProposed) ;
            likelihoodRatio = exp(lLikelihoodProposed - lLikelihood) ;
            
            
            acceptance = min(1, likelihoodRatio * priorRatio * proposalRatio);
            accRatio = accRatio + acceptance/iterations ;
            if (rand(1) <= acceptance) %&& swDistAccept
                lLikelihood = lLikelihoodProposed;
                lLikelihoods = lLikelihoodsProposed;
                regProfileActivationAll = regProfileActivationAllProposed;
                regProfileThreshold = regProfileThresholdProposed;
                if (strcmp(thresholdPrior{1},'gradient'))
                    priorThresholds = priorThresholdsProposed ;
                end
                
                accept(1) = accept(1) + 1 ;
            end
        end
            
        
        
        %%% do the S (swap) move
    elseif(moveIdx == 2)

        %%% add another S jump
        proposalNum(2) = proposalNum(2) + 1 ;
        
        %%% pick one of the selected regs
        remSelParentIdx = randi(numSelRegs);
        remParentIdx = selRegs(remSelParentIdx) ;
        
        %%% pick one of the non-selected regs
        addNonSelParentIdx = randi(numNonSelRegs);
        addParentIdx = nonSelRegs(addNonSelParentIdx) ;
        
        %%% swap them
        selRegsProposed = selRegs(selRegs ~= remParentIdx) ;
        selRegsProposed = sort([selRegsProposed addParentIdx]) ;
        numSelRegsProposed = length(selRegsProposed) ;
        
        % get the current thresholds
        regProfileThresholdProposed = regProfileThreshold;
        
        % sample a new threshold for the new reg
        % note that we only use the threshold for the selected reg so no
        % need to pick the thresholds of selected regs
        if (strcmp(proposalsType,'normal'))
            th = fastnormrnd(regProfileMean(addParentIdx),rwSd,1) ;
            while (th <= minRegProfiles(addParentIdx) || th >= maxRegProfiles(addParentIdx))
                th = fastnormrnd(regProfileMean(addParentIdx),rwSd,1) ;
            end
        elseif (strcmp(proposalsType,'uniform'))
            th = minRegProfiles(addParentIdx) + rand*(maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx)) ;
        end
        regProfileThresholdProposed(addParentIdx) = th ;
        
        priorRatioNumSelRegs = 1 ;
        if (strcmp(regPrior{1},'range') )
            priors = regPrior{2} ;
            if (numSelRegsProposed>0)
                [~,idxCombiSelRegsProposed,~] = intersect(nchoosek(allRegs,numSelRegsProposed),selRegsProposed,'rows') ;
                priorSelRegsProposed = priors{numSelRegsProposed}(idxCombiSelRegsProposed,1) ;
            else
                priorSelRegsProposed = 1 ;
            end
            priorRatioSelRegs = priorSelRegsProposed/priorSelRegs ;
        else
            priorSelRegsProposed = 0 ;
            priorRatioSelRegs = 1 ;
        end
        
        if (strcmp(thresholdPrior{1},'gradient') )
            priorThresholdsProposed = priorThresholds ;
            idxTh = find(th >= levels(:,addParentIdx),1,'last' ) ;
            priorThresholdsProposed(addParentIdx) = thresholdPrior{2}(idxTh,addParentIdx) ;
            priorRatioThresholds = priorThresholdsProposed(addParentIdx)/priorThresholds(remParentIdx) ;
        else
            priorRatioThresholds = (maxRegProfiles(remParentIdx) - minRegProfiles(remParentIdx))/(maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx)) ;
        end
        
        priorRatio = priorRatioNumSelRegs*priorRatioSelRegs*priorRatioThresholds ;
        
        if (priorRatio ~= 0 && ~isnan(priorRatio) && ~isinf(priorRatio))
            
            for x = 1:numExps
                
                time = timeAll{x};
                
                % get the activation for the selected regs
                regProfileActivationSelRegsProposed = regProfileActivationAll{x}(:,selRegsProposed);
                
                % change the activation of the new reg
                % first get its expression
                regProfileAddParent = regProfileAll{x}(:,addParentIdx);
                
                % then get the activation
                regProfileActivationSelRegsProposed(hSteps(x):end,selRegsProposed==addParentIdx) = getRegProfileActivation(regProfileAddParent(hSteps(x):end),th);
                % set the act func before delay as the first after delay to avoid switches
                regProfileActivationSelRegsProposed(1:hSteps(x)-1,selRegsProposed==addParentIdx) = regProfileActivationSelRegsProposed(hSteps(x),selRegsProposed==addParentIdx) ;
                
                regProfileActivationAllProposed{x}(:,addParentIdx) = regProfileActivationSelRegsProposed(:,selRegsProposed==addParentIdx);
                
                twosMat = repmat(pow2(0:(numSelRegs-1)),timeFuncSize,1);
                stateChanges = getGeneActivation(regProfileActivationSelRegsProposed,twosMat);
                [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x));
                switchesAll{x} = switches;
                statesAll{x} = states;
            end
            
            statesUsed = getStatesUsed(statesAll) ;
            
            [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,replicatesAll,lsquares);
            m0All(m0All<0)=0;
            birthRates(birthRates<0)=0;

                for x = 1:numExps
                    
                    br = zeros(length(statesAll{x}),1) ;
                    for y=1:length(br)
                        br(y) = find(statesAll{x}(y) == statesUsed) ;
                    end
                    targetmRNAProposed = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                    
                    %calculate SSE of profile to data
                    e_hat{x} = calculateEHat(targetmRNAProposed,childDataAll{x},eHatMatAll{x},replicatesAll(x));
                    
                    if strcmp(lsquares.type,'wls')
                        lLikelihoodsProposed(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
                    else
                        lLikelihoodsProposed(x) = logLikelihood(e_hat{x},prec(x))  ;
                    end
                end                
                lLikelihoodProposed = sum(lLikelihoodsProposed) ;
                likelihoodRatio = exp(lLikelihoodProposed - lLikelihood) ;
                
                if (strcmp(proposalsType,'normal'))
                    proposalRatio = normpdf(regProfileThreshold(remParentIdx),regProfileMean(remParentIdx),rwSd)/normpdf(regProfileThresholdProposed(addParentIdx),regProfileMean(addParentIdx),rwSd) ;
                elseif (strcmp(proposalsType,'uniform'))
                    proposalRatio = (maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx))/(maxRegProfiles(remParentIdx) - minRegProfiles(remParentIdx)) ;
                end
                
                acceptance = min(1, likelihoodRatio * priorRatio * proposalRatio);
                accRatio = accRatio + acceptance/iterations ;
                if (rand(1) <= acceptance) 
                    lLikelihood = lLikelihoodProposed;
                    lLikelihoods = lLikelihoodsProposed;
                    regProfileActivationAll = regProfileActivationAllProposed;
                    regProfileThreshold = regProfileThresholdProposed;
                    
                    selRegs = selRegsProposed ;
                    nonSelRegs = setdiff(allRegs,selRegs) ;
                    
                    priorSelRegs = priorSelRegsProposed ;
                    if (strcmp(thresholdPrior{1},'gradient'))
                        priorThresholds = priorThresholdsProposed ;
                    end
                    
                    accept(2) = accept(2) + 1 ;
                end
            
        end
        
        
        %%% do the add reg (B) move
    elseif(moveIdx == 3)

        %%% add another B jump
        proposalNum(3) = proposalNum(3) + 1 ;
        
        
        %%% pick one of the non-selected regs
        addNonSelParentIdx = randi(numNonSelRegs);
        addParentIdx = nonSelRegs(addNonSelParentIdx) ;
        
        %%% add it to selected Regs
        selRegsProposed = sort([selRegs addParentIdx]) ;
        numSelRegsProposed = numSelRegs + 1 ;
        
        %%% get the thresholds
        % call the current
        regProfileThresholdProposed = regProfileThreshold;
        
        % change the one we added using truncated normal
        if (strcmp(proposalsType,'normal'))
            th = fastnormrnd(regProfileMean(addParentIdx),rwSd,1) ;
            while (th <= minRegProfiles(addParentIdx) || th >= maxRegProfiles(addParentIdx))
                th = fastnormrnd(regProfileMean(addParentIdx),rwSd,1) ;
            end
        elseif (strcmp(proposalsType,'uniform'))
            th = minRegProfiles(addParentIdx) + rand*(maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx)) ;
        end
        regProfileThresholdProposed(addParentIdx) = th ;
        
        %%% this uses a prior vaguely centered in the middle of the
        %%% expression levels.
        priorRatioNumSelRegs = numRegsPrior(numSelRegsProposed+1)/numRegsPrior(numSelRegs+1) ; 
        if strcmp(regPrior{1},'range')
            priors = regPrior{2} ;
            if (numSelRegsProposed>0)
                [~,idxCombiSelRegsProposed,~] = intersect(nchoosek(allRegs,numSelRegsProposed),selRegsProposed,'rows') ;
                priorSelRegsProposed = priors{numSelRegsProposed}(idxCombiSelRegsProposed,1) ;
            else
                priorSelRegsProposed = 1 ;
            end
            priorRatioSelRegs = priorSelRegsProposed/priorSelRegs ;
        else
            priorSelRegsProposed = 0 ;
            % the above is not the correct prior but simply to get an input to go into update in case of uniform prior
            priorRatioSelRegs = numSelRegsProposed/(numRegs - numSelRegs) ;
        end
        
        if (strcmp(thresholdPrior{1},'gradient'))
            priorThresholdsProposed = priorThresholds ;
            idxTh = find(th >= levels(:,addParentIdx),1,'last' ) ;
            priorThresholdsProposed(addParentIdx) = thresholdPrior{2}(idxTh,addParentIdx) ;
            priorRatioThresholds = priorThresholdsProposed(addParentIdx) ;
        else
            priorRatioThresholds = 1/(maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx)) ;
        end
        
        proposalMoveChoiceRatio = truncPoisson(numSelRegs,numSelRegsPropLambda,maxNumRegs)/truncPoisson(numSelRegsProposed,numSelRegsPropLambda,maxNumRegs) ;
        proposalParentChoiceRatio = (numRegs - numSelRegs)/numSelRegsProposed ;
        if (strcmp(proposalsType,'normal'))
            proposalThresholdMoveRatio = 1/normpdf(regProfileThresholdProposed(addParentIdx),regProfileMean(addParentIdx),rwSd) ;
        elseif (strcmp(proposalsType,'uniform'))
            proposalThresholdMoveRatio = (maxRegProfiles(addParentIdx) - minRegProfiles(addParentIdx)) ;
        end
        proposalRatio = proposalMoveChoiceRatio*proposalParentChoiceRatio*proposalThresholdMoveRatio ;
        priorRatio = priorRatioNumSelRegs*priorRatioSelRegs*priorRatioThresholds   ;
        
        if (priorRatio ~= 0 && ~isnan(priorRatio) && ~isinf(priorRatio))
            
            for x = 1:numExps
                
                time = timeAll{x};
                
                % get the current activation function for the selected regs
                regProfileActivationSelRegsProposed = regProfileActivationAll{x}(:,selRegsProposed);
                
                % change the activation function of the added reg
                % first get its expression
                regProfileAddParent = regProfileAll{x}(:,addParentIdx);
                
                regProfileActivationSelRegsProposed(hSteps(x):end,selRegsProposed==addParentIdx) = getRegProfileActivation(regProfileAddParent(hSteps(x):end),th);
                % set the act func before delay as the first after delay to avoid switches
                regProfileActivationSelRegsProposed(1:hSteps(x)-1,selRegsProposed==addParentIdx) = regProfileActivationSelRegsProposed(hSteps(x),selRegsProposed==addParentIdx) ;
                
                regProfileActivationAllProposed{x}(:,addParentIdx) = regProfileActivationSelRegsProposed(:,selRegsProposed==addParentIdx);
                
                twosMat = repmat(pow2(0:(numSelRegsProposed-1)),timeFuncSize,1);
                stateChanges = getGeneActivation(regProfileActivationSelRegsProposed,twosMat);
                [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x));
                switchesAll{x} = switches;
                statesAll{x} = states;
            end
            
            statesUsed = getStatesUsed(statesAll);
            
            [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,replicatesAll,lsquares);
            m0All(m0All<0)=0;
            birthRates(birthRates<0)=0;

            for x = 1:numExps
                
                br = zeros(length(statesAll{x}),1) ;
                for y=1:length(br)
                    br(y) = find(statesAll{x}(y) == statesUsed) ;
                end
                targetmRNAProposed = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                
                %calculate SSE of profile to data
                e_hat{x} = calculateEHat(targetmRNAProposed,childDataAll{x},eHatMatAll{x},replicatesAll(x)) ;
                
                if strcmp(lsquares.type,'wls')
                    lLikelihoodsProposed(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
                else
                    lLikelihoodsProposed(x) = logLikelihood(e_hat{x},prec(x))  ;
                end
                
            end
            lLikelihoodProposed = sum(lLikelihoodsProposed) ;
            likelihoodRatio = exp(lLikelihoodProposed - lLikelihood) ;
            acceptance = min(1, likelihoodRatio * priorRatio * proposalRatio);
            accRatio = accRatio + acceptance/iterations ;
            if (rand(1) <= acceptance)
                lLikelihood = lLikelihoodProposed;
                lLikelihoods = lLikelihoodsProposed;
                regProfileActivationAll = regProfileActivationAllProposed;
                regProfileThreshold = regProfileThresholdProposed;
                
                selRegs = selRegsProposed ;
                numSelRegs = numSelRegsProposed ;
                nonSelRegs = setdiff(allRegs,selRegs) ;
                numNonSelRegs = numRegs - numSelRegs ;
                
                priorSelRegs = priorSelRegsProposed ;
                
                if (strcmp(thresholdPrior{1},'gradient'))
                    priorThresholds = priorThresholdsProposed ;
                end
                
                accept(3) = accept(3) + 1 ;
            end
            
        end
        
        %%% do the remove reg (D) move
    else
        %%% add another D jump
        proposalNum(4) = proposalNum(4) + 1 ;
        
        %%% pick one of the selected regs
        remSelParentIdx = randi(numSelRegs);
        remParentIdx = selRegs(remSelParentIdx) ;
        
        selRegsProposed = selRegs(selRegs ~= remParentIdx) ;
        numSelRegsProposed = numSelRegs - 1 ;
                
        priorRatioNumSelRegs = numRegsPrior(numSelRegsProposed+1)/numRegsPrior(numSelRegs+1) ; % truncPoisson(numSelRegsProposed,numSelRegsPriorLambda,maxNumRegs)/truncPoisson(numSelRegs,numSelRegsPriorLambda,maxNumRegs)
        if (strcmp(regPrior{1},'range') )
            priors = regPrior{2} ;
            if (numSelRegsProposed>0)
                [~,idxCombiSelRegsProposed,~] = intersect(nchoosek(allRegs,numSelRegsProposed),selRegsProposed,'rows') ;
                priorSelRegsProposed = priors{numSelRegsProposed}(idxCombiSelRegsProposed,1) ;
            else
                priorSelRegsProposed = 1 ;
            end
            priorRatioSelRegs = priorSelRegsProposed/priorSelRegs ;
        else
            priorSelRegsProposed = 0 ;
            priorRatioSelRegs = (numRegs - numSelRegsProposed)/numSelRegs ;
        end
        
        if (strcmp(thresholdPrior{1},'gradient'))
            priorRatioThresholds = 1/priorThresholds(remParentIdx) ;
        else
            priorRatioThresholds = (maxRegProfiles(remParentIdx) - minRegProfiles(remParentIdx)) ;
        end
        priorRatio = priorRatioNumSelRegs*priorRatioSelRegs*priorRatioThresholds ;
        
        if (priorRatio ~= 0 && ~isnan(priorRatio) && ~isinf(priorRatio))
            
            
            for x = 1:numExps
                
                time = timeAll{x};
                
                % get the activation function of the selected regs
                regProfileActivationSelRegsProposed = regProfileActivationAll{x}(:,selRegsProposed);
                
                % no need to change anything in act function. I only removed
                % one column
               
                twosMat = repmat(pow2(0:(numSelRegsProposed-1)),timeFuncSize,1);
                stateChanges = getGeneActivation(regProfileActivationSelRegsProposed,twosMat);
                [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x));
                switchesAll{x} = switches;
                statesAll{x} = states;
            end
            
            statesUsed = getStatesUsed(statesAll);
            
            [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,replicatesAll,lsquares);
            m0All(m0All<0)=0;
            birthRates(birthRates<0)=0;

            for x = 1:numExps
                
                br = zeros(length(statesAll{x}),1) ;
                for y=1:length(br)
                    br(y) = find(statesAll{x}(y) == statesUsed) ;
                end
                targetmRNAProposed = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                
                %calculate SSE of profile to data
                e_hat{x} = calculateEHat(targetmRNAProposed,childDataAll{x},eHatMatAll{x},replicatesAll(x));
                if strcmp(lsquares.type,'wls')
                    lLikelihoodsProposed(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
                else
                    lLikelihoodsProposed(x) = logLikelihood(e_hat{x},prec(x))  ;
                end
            end
            lLikelihoodProposed = sum(lLikelihoodsProposed) ;
            likelihoodRatio = exp(lLikelihoodProposed - lLikelihood) ;
            
            proposalMoveChoiceRatio = truncPoisson(numSelRegs,numSelRegsPropLambda,maxNumRegs)/truncPoisson(numSelRegsProposed,numSelRegsPropLambda,maxNumRegs) ;
            proposalParentChoiceRatio = numSelRegs/(numRegs - numSelRegsProposed) ;
            if (strcmp(proposalsType,'normal'))
                proposalThresholdMoveRatio = normpdf(regProfileThresholdProposed(remParentIdx),regProfileMean(remParentIdx),rwSd) ;
            elseif (strcmp(proposalsType,'uniform'))
                proposalThresholdMoveRatio = 1/(maxRegProfiles(remParentIdx) - minRegProfiles(remParentIdx)) ;
            end
            
            proposalRatio = proposalMoveChoiceRatio*proposalParentChoiceRatio*proposalThresholdMoveRatio ;
            
            acceptance = min(1, likelihoodRatio * priorRatio * proposalRatio);
            accRatio = accRatio + acceptance/iterations ;
            if (rand(1) <= acceptance) %&& swDistAccept
                lLikelihood = lLikelihoodProposed;
                lLikelihoods = lLikelihoodsProposed;
                regProfileActivationAll = regProfileActivationAllProposed;
                regProfileThreshold = regProfileThresholdProposed;
                
                selRegs = selRegsProposed ;
                numSelRegs = numSelRegsProposed ;
                nonSelRegs = setdiff(allRegs,selRegs) ;
                numNonSelRegs = numRegs - numSelRegs ;
                
                priorSelRegs = priorSelRegsProposed ;
                
                accept(4) = accept(4) + 1 ;
            end
            
        end
    end
        
    if deg_sampled
        %%% lets change the mRNA degradation rate next, as the other parameters are dependent on this
        degProposed = fastnormrnd(deg,deg_rwSd,1) ; % exp(log(deg) + fastnormrnd(0,0.25,1));
        while (degProposed < minDegRate || maxDegRate < degProposed)
            degProposed = fastnormrnd(deg,0.25,1);
        end
        for x = 1:numExps
            
            regProfileActivationSelRegs = regProfileActivationAll{x}(:,selRegs);
            time = timeAll{x};
            twosMat = repmat(pow2(0:(numSelRegs-1)),timeFuncSize,1) ;
            stateChanges = getGeneActivation(regProfileActivationSelRegs,twosMat) ;
            [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x)) ;
            switchesAll{x} = switches;
            statesAll{x} = states;
        end
        statesUsed = getStatesUsed(statesAll);
        
        [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,degProposed,YAll,ntimepointsAll,replicatesAll,lsquares);
        m0All(m0All<0)=0;
        birthRates(birthRates<0)=0;
        
        for x = 1:numExps
            br = zeros(length(statesAll{x}),1) ;
            for y=1:length(br)
                br(y) = find(statesAll{x}(y) == statesUsed) ;
            end
            targetmRNAProposed = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),degProposed,switchesAll{x});
            
            e_hat{x} = calculateEHat(targetmRNAProposed,childDataAll{x},eHatMatAll{x},replicatesAll(x));
            if strcmp(lsquares.type,'wls')
                lLikelihoodsProposed(x) = logLikelihoodWLS(e_hat{x},prec(x),lsquares.Q{x})  ;
            else
                lLikelihoodsProposed(x) = logLikelihood(e_hat{x},prec(x))  ;
            end
        end
        lLikelihoodProposed = sum(lLikelihoodsProposed) ;
        likelihoodRatio = exp(lLikelihoodProposed - lLikelihood) ;
        
        priorRatio = gampdf(degProposed,prior_deg_a0,prior_deg_b0)/gampdf(deg,prior_deg_a0,prior_deg_b0);
        proposalRatio = 1 ;
        
        acceptance = min(1, likelihoodRatio * priorRatio * proposalRatio);
        if rand(1) <= acceptance
            lLikelihood = lLikelihoodProposed;
            lLikelihoods = lLikelihoodsProposed;
            deg = degProposed;
            accept(5) = accept(5) + 1;
        end
        
    end
    
    lLikelihoodHistory(i,1:numExps) = lLikelihoods';
    lLikelihoodHistory(i,end) = lLikelihood;
    likelihoodRatioHistory(i) = likelihoodRatio ;
    
    degHistory(i) = deg ;
    regProfileThresholdHistory(i,selRegs) = regProfileThreshold(selRegs)';
    regProfileThresholdHistory(i,nonSelRegs) = nan ;
    
    m0History(i,:) = m0All';
    precHistory(i,:) = prec';
    selRegsHistory{i} = selRegs ;
    numSelRegsHistory(i) = numSelRegs ;
    
    for z=1:numRegs
        selParentHistoryPc(z) = selParentHistoryPc(z) + sum(selRegs == z)/iterations ;
        numSelRegsHistoryPc(z) = numSelRegsHistoryPc(z) + sum(numSelRegs == z)/iterations ;
    end
    
    noFail = sum(birthRates < 0) == 0 ;
    if ~noFail
        keyboard
    end
    moveIdxHistory(i) = moveIdx ;
    birthratesHistory{i} = birthRates ;
    for x = 1:numExps
        switchesHistory{i,x} = switchesAll{x} ;
        statesHistory{i,x} = statesAll{x} ;
    end    
    if (strcmp(lsquares.type,'wls'))
        lsquaresWeightsPhiHistory(i,:) = lsquares.weightsPhi ;            
    end
end
disp(['acceptance ratio=' num2str(accRatio)])
disp(['acceptance ratio for move M=' num2str(accept(1)/proposalNum(1))])
disp(['acceptance ratio for move S=' num2str(accept(2)/proposalNum(2))])
disp(['acceptance ratio for move B=' num2str(accept(3)/proposalNum(3))])
disp(['acceptance ratio for move D=' num2str(accept(4)/proposalNum(4))])
disp(['acceptance ratio for move deg=' num2str(accept(5)/iterations)])

boxOn = 0 ;
if mcmcPlotsOn
    figure('Name','Parameter chains','NumberTitle','off','WindowStyle','docked');
    subplot('Position',[0.1 0.6 0.35 0.35]);
    plot(degHistory);
    title('degradation rate','interpreter','none','FontSize',12.5,'FontWeight','bold','linewidth',2);
    set(gca,'XTickMode','auto', 'FontSize',12,'FontWeight','bold')
    axis tight
    
    subplot('Position',[0.6 0.6 0.35 0.35]);
    leg = cell(numExps,1) ;
    for x = 1:numExps
        plot(precHistory(:,x),'linewidth',2);
        hold all;
        set(gca,'XTickMode','auto', 'FontSize',12,'FontWeight','bold')
        leg{x} = char(['exp' num2str(x)]) ;
    end
    ylim([0 max(prctile(precHistory(100:end,:),99.5))+ max(prctile(precHistory(100:end,:),0.05))])
    xlim([0 iterations])
    legend(leg,'location','best')
    title('precision','FontSize',12.5,'FontWeight','bold');
    s=subplot('Position',[0.1 0.1 0.6 0.35]);
    leg = cell(numRegs,1) ;
    for z = 1:numRegs
        pl(z) = plot(regProfileThresholdHistory(:,z),'linewidth',2);
        hold all
        leg{z} = char(['TF' num2str(z)]) ;
    end
    l = legend(leg,'location','northwest');
    title('TF-threshold level','FontSize',12.5,'FontWeight','bold');
    set(gca,'XTickMode','auto', 'FontSize',12,'FontWeight','bold')
    hold off;
    axis tight
    if boxOn
        p = get(s,'Position') ;
        a = axes('position',[p(1)+p(3)-0.125 p(2)+0.05 .1 .1]);
        set(a,'box','on')
        [~,zidx] = min(sum(isnan(regProfileThresholdHistory))) ;
        nans=1 ;
        k=1 ;
        while nans
            idxs = iterations-50*k:iterations-(k-1)*50 ;
            if sum(isnan(regProfileThresholdHistory(idxs,zidx)))==0
                nans=0 ;
            else
                k=k+1;
            end
        end
        col = get(pl(zidx),'Color') ;
        plot(idxs,regProfileThresholdHistory(idxs,zidx),'Color',col,'linewidth',2)
        axis tight
        set(a,'FontSize',10,'FontWeight','bold')
        ytick = get(a,'YTick') ;
        set(a,'YTick',[ytick(1) ytick(end)] )
        title(['TF' num2str(zidx) ' magnified'])
    end
    s=subplot('Position',[0.8 0.1 0.15 0.35]);
    boxplot(lLikelihoodHistory(:,numExps+1));
    set(gca,'FontSize',12.5,'FontWeight','bold')
    set(gca,'XTick',[])
    title('log likelihood','FontSize',12.5,'FontWeight','bold');
    YTick = get(gca,'YTick') ;
    if length(YTick)>4
        YTick = YTick(1:2:end) ;
    end
    set(gca,'YTick',YTick) ;
end


%% outputs
output = struct;
output.allRegs = allRegs ;
output.numSelRegsHistory =  numSelRegsHistory ;
output.selRegsHistory = selRegsHistory ;
output.regProfileThresholdHistory = regProfileThresholdHistory ;
output.switchesHistory = switchesHistory ;
output.statesHistory = statesHistory ;
output.birthratesHistory = birthratesHistory ;
output.m0History = m0History;

output.degHistory = degHistory ;
output.precHistory = precHistory ;
output.lLikelihoodHistory = lLikelihoodHistory ;
output.likelihoodRatioHistory = likelihoodRatioHistory ;

output.moveIdxHistory = moveIdxHistory ;
output.accept = accept ;
output.proposalNum = proposalNum ;
output.accRatio = accRatio ;

output.numSelRegsHistoryPc = numSelRegsHistoryPc ;
output.selParentHistoryPc = selParentHistoryPc ;

output.lsquares = lsquares ;
output.lsquaresWeightsPhiHistory = lsquaresWeightsPhiHistory ;

warning(origWarnings);

end




function [] = boxplots(data)
%each column are samples for box plot
xticklabel = cell(size(data,2),1) ;
for x = 1:size(data,2)
    d = sort(data(:,x),'ascend');
    p = getPercentiles([2.5 25 50 75 97.5],d);
    
    plot([x*2-0.25 x*2+0.25],[p(1) p(1)],'-k');
    hold on;
    plot([x*2 x*2],[p(1) p(2)],'--k');
    plot([x*2-0.5 x*2+0.5],[p(2) p(2)],'-k');
    plot([x*2-0.5 x*2+0.5],[p(3) p(3)],'-k');
    plot([x*2-0.5 x*2+0.5],[p(4) p(4)],'-k');
    plot([x*2-0.5 x*2-0.5],[p(2) p(4)],'-k');
    plot([x*2+0.5 x*2+0.5],[p(2) p(4)],'-k');
    plot([x*2 x*2],[p(4) p(5)],'--k');
    plot([x*2-0.25 x*2+0.25],[p(5) p(5)],'-k');
    xticklabel{x} = num2str(x) ;
end
xlim([1 size(data,2)*2+1]);
set(gca,'XTick',2:2:2*size(data,2))
set(gca,'XTickLabel',xticklabel)
hold off;
end

function [val] = getPercentiles(percentiles,data)
n = length(data);
val = data(floor((n / 100) .* percentiles + 0.5));
end

function prob = truncPoisson(k,lambda,kmax)

if (k > -1 && k < kmax+1)
    prob = exp(-lambda)*lambda^k/factorial(k) ;
else
    prob = 0 ;
end
%%% note that we do not compute the normalising constant for the truncation
%%% we don't bcs it's not necessary in the rates but be careful in case you
end

function [states] = getGeneActivation(regProfile,twos)

states = binaryMatToRow(regProfile,twos);

end

function [rows] = binaryMatToRow(binMatrix,twos)

rows = sum(twos .* binMatrix,2)+1;

end

function [activation] = getRegProfileActivation(expression,threshold)

if (isempty(expression) == 0)
    activation = expression >= threshold;
else
    [l,~] = size(expression) ;
    activation = zeros(l,1) ;
end

end

function [e_hat] = calculateEHat(profile,data,eHatMat,nReps)
for x = 1:nReps
    eHatMat(x,:) = profile ; % create the fit matrix with nrows = replicates, ncolumns = timepoints
end
e_hat = eHatMat - data ;% compute the residuals
e_hatT = e_hat' ;
e_hat = e_hatT(:) ;
end

function r = fastnormrnd(mu,sigma,numRands)
r = randn(1,numRands) .* sigma + mu;
end


function statesUsed = getStatesUsed(states)

statesUsed = unique(states{1});

for x = 2:length(states)
    statesUsed = unique([states{x} statesUsed]);
end

end

function llik = logLikelihood(e_hat,prec)

sig2 = 1 / prec ;
llikvec = log(sqrt(2*pi*sig2)) + ((e_hat).^2)./(2*sig2);
llik = -1 * sum(llikvec);

end

function llik = logLikelihoodWLS(e_hat,prec,Q)

sig2 = 1 / prec ;
Sig = sig2*Q ;
n = length(e_hat) ;
llik = - n*log(2*pi)/2 - n*log(sig2)/2 - sum(log(diag(Q)))/2 - (e_hat'*(Sig\e_hat))/2 ;

end









