function [output] =  createFinalOutput(profilesOutput,priorsPropOutput,mcmcOutput,varargin)

% output
% freq of numRegs 
% freq of regs given numRegs
% freq of regs
% freq of regs(1), freq of regs(2)
% histogram of thresholds for numRegs,regs
% activation function and fit for mode of numRegs,regs

%% default variables
regsFreq = true ;
interactivePlots = 0 ;
finalPlotsOn = 1 ;
regNames = [];
plot2modelsOn = 0 ;
ratesFreqOn = 1 ;

%% process optional inputs
if mod(length(varargin),2)
    error('Syntax error on command line. Enter Name Value pairs');
end
for i=1:2:length(varargin)-1
    switch varargin{i}        
        case ('regulatorsFreq')
            regsFreq = varargin{i+1};
        case {'finalPlotsOn'}
            finalPlotsOn = varargin{i+1};       
        case ('interactivePlots')
            interactivePlots = varargin{i+1}; 
        case('regulatorNames')
            regNames = varargin{i+1};
        case('plot2modelsOn')
            plot2modelsOn = varargin{i+1};
        otherwise
            error(['Unknown property ' varargin{i}  ' in the command line']);
    end
end

%% get various inputs
experiments = profilesOutput.experiments ;

selRegsHistory = mcmcOutput.selRegsHistory ; % cell with entries the vectors of selRegs in each iteration
numSelRegsHistory = mcmcOutput.numSelRegsHistory ; % vector with entries the number of selRegs in each iteration

iterations = length(selRegsHistory) ;

numRegs = profilesOutput.nregulators ; % max num Regs allowed in MCMC
maxNumRegs = priorsPropOutput.maxNumRegs ; % max num Regs allowed in MCMC
numExps = profilesOutput.nexps;
timeFuncSize = profilesOutput.timeFuncSize ;
minRegProfiles = profilesOutput.minRegProfiles ;
maxRegProfiles = profilesOutput.maxRegProfiles ;
hSteps = profilesOutput.hSteps ;
minSwTime = profilesOutput.minSwTime  ;
eHatMatAll = profilesOutput.eHatMatAll ;

regProfileAll = profilesOutput.regProfileAll ;
degHistory = mcmcOutput.degHistory ;
precHistory = mcmcOutput.precHistory ;
regProfileThresholdHistory = mcmcOutput.regProfileThresholdHistory ;


if ratesFreqOn
    m0History = mcmcOutput.m0History ;
    birthRatesHistory = mcmcOutput.birthratesHistory ;
    switchesHistory = mcmcOutput.switchesHistory ;
end
lsquaresWeightsPhiHistory = mcmcOutput.lsquaresWeightsPhiHistory ;

minDegRate = priorsPropOutput.priors.minDegRate ;
maxDegRate = priorsPropOutput.priors.maxDegRate ;

targetDataAll = profilesOutput.targetDataAll ;

timepointsAll = profilesOutput.timepointsAll ;
ntimepointsAll = profilesOutput.ntimepointsAll ;
nreplicatesAll = profilesOutput.nreplicatesAll ;

minTargetData = profilesOutput.minTargetData ;
maxTargetData = profilesOutput.maxTargetData ;

profileFuncSize = priorsPropOutput.profileFuncSize ;

if isempty(regNames)
    regNames = profilesOutput.regulatorNames ;
end

% how many regs?
% freq of number of regs
numRegsFreq = mcmcOutput.numSelRegsHistoryPc ;
numRegsFreq(end+1:numRegs) = 0;

output = struct ;
output.numRegsFreq = numRegsFreq ;

if regsFreq
    % which reg? pair of regs? triplet of regs? etc.
    % First given the number of regs included
    regsCombiGivenNumRegsFreq = cell(maxNumRegs,1) ; % cell with entries the freq and (combi) of regs
    for q = 1:min(maxNumRegs,numRegs)
        leCombis =  nchoosek(numRegs,q) ; % how many combis we have
        regsCombiGivenNumRegsFreq{q} = zeros(leCombis,q+1) ;  % we need that many freq and also columns including the combi
        regsCombiGivenNumRegsFreqQ = regsCombiGivenNumRegsFreq{q} ;
        regsCombiGivenNumRegsFreqQ(:,2:end) = nchoosek(1:numRegs,q) ; % insert the combis
        idxNumSelRegsQ = find(numSelRegsHistory == q) ; % find which iterations had numRegs=q
        leNumSelRegsQ = length(idxNumSelRegsQ) ;  % how many they are
        if (leNumSelRegsQ ~= 0) % if we have at least one iter with q regs
            selRegsHistoryGivenNumRegs = zeros(leNumSelRegsQ,q) ;  % create a matrix with elements the regs for all iterations with numRegs=q
            for i = 1:leNumSelRegsQ
                selRegsHistoryGivenNumRegs(i,:) = selRegsHistory{idxNumSelRegsQ(i)}  ;
            end
            for j = 1:leCombis
                combi = regsCombiGivenNumRegsFreqQ(j,2:end); % get the combi
                % find how many rows are equal to combi, divide with all iterations
                regsCombiGivenNumRegsFreqQ(j,1) = length(find(ismember(selRegsHistoryGivenNumRegs,combi,'rows')==1))/leNumSelRegsQ ;
            end
        end
        regsCombiGivenNumRegsFreq{q} = regsCombiGivenNumRegsFreqQ ;
    end
    
    % secondly, marginals: which single/pair of/triplet of reg(s) is picked most across all chosen models?
    regsCombiFreq = cell(maxNumRegs,1) ; 
    for q = 1 %:maxNumRegs
        leCombis =  nchoosek(numRegs,q) ; 
        regsCombiFreq{q} = zeros(leCombis,q+1) ; 
        regsCombiFreqQ = regsCombiFreq{q} ;
        regsCombiFreqQ(:,2:end) = nchoosek(1:numRegs,q) ;
        for j=1:leCombis
            combi = regsCombiGivenNumRegsFreq{q}(j,2:end) ;% get the combi
            regsCombiFreqQ(j,1) = regsCombiGivenNumRegsFreq{q}(j,1)*numRegsFreq(q) ; % first get the freq when numPar = length(combi)
            for g = q+1:min(maxNumRegs,numRegs) % for larger numPar
                combis = regsCombiGivenNumRegsFreq{g}(:,2:end) ; % get all the combis
                regsCombiFreqQ(j,1) = regsCombiFreqQ(j,1) + sum(regsCombiGivenNumRegsFreq{g}(sum(ismember(combis,combi),2) == q,1)*numRegsFreq(g));
            end
            % the final sum has taken the freq*numRegsFreq of all combis
            % having combi as subset
        end
        regsCombiFreq{q} = regsCombiFreqQ ;
    end
    
    % which model? (joint of number and combi of regs)
    % the freq of each model within the whole chain
    %%% there must be a faster way to do this but I used the initial code
    %%% you only need the freq of all observed. Ignoring the numbers etc.
    modelFreq = cell(maxNumRegs,1) ;
    for q = 1:min(maxNumRegs,numRegs)
        regsCombiGivenNumRegsFreqQ = regsCombiGivenNumRegsFreq{q} ;
        if ~isempty(regsCombiGivenNumRegsFreqQ)
            modelFreqQ = regsCombiGivenNumRegsFreqQ ;
            modelFreqQ(:,1) = modelFreqQ(:,1)*sum(mcmcOutput.numSelRegsHistory==q)/iterations ;
            [~,i] = sort(modelFreqQ(:,2));
            modelFreq{q} = modelFreqQ(i,:) ;
            
        end
    end
    
    % thresholds and deg for given model = (numRegs,regs)
    thresholdsGivenModel = cell(maxNumRegs,1) ;
    thresholdsGivenModelKsDens = cell(maxNumRegs,1) ;
    degGivenModel = cell(maxNumRegs,1) ;
    degGivenModelKsDens = cell(maxNumRegs,1) ;
    precGivenModel = cell(maxNumRegs,1) ;
    precGivenModelKsDens = cell(maxNumRegs,1) ;
    lsquaresWeightsPhiGivenModel = cell(maxNumRegs,1) ;
    if ratesFreqOn
        birthratesGivenModel = cell(maxNumRegs,1) ;
        m0GivenModel = cell(maxNumRegs,1) ;
        switchesGivenModel = cell(maxNumRegs,1) ;
    end
    for q = 1:min(maxNumRegs,numRegs)  
        leCombis = size(regsCombiGivenNumRegsFreq{q},1) ;
        
        thresholdsGivenModel{q} = cell(leCombis,2) ;
        thresholdsGivenModelKsDens{q} = cell(leCombis,2) ;
        degGivenModel{q} = cell(leCombis,2) ;
        degGivenModelKsDens{q} = cell(leCombis,2) ;
        precGivenModel{q} = cell(leCombis,2) ;
        precGivenModelKsDens{q} = cell(leCombis,2) ;
        lsquaresWeightsPhiGivenModel{q} = cell(leCombis,2) ;
        if ratesFreqOn
            birthratesGivenModel{q} = cell(leCombis,2) ;
            m0GivenModel{q} = cell(leCombis,2) ;
            switchesGivenModel{q} = cell(leCombis,2) ;
        end
        idxNumSelRegsQ = find(numSelRegsHistory == q) ; % find which iterations had numRegs=q
        leNumSelRegsQ = length(idxNumSelRegsQ) ;  % how many they are
        selRegsHistoryGivenNumRegs = zeros(leNumSelRegsQ,q) ;  % create a matrix with elements the regs for all iterations with numRegs=q
        for i = 1:leNumSelRegsQ
            selRegsHistoryGivenNumRegs(i,:) = selRegsHistory{idxNumSelRegsQ(i)}  ;
        end
        idxNumRegCombiNotEmpty = zeros(leCombis,1) ;
        for j=1:leCombis
            combi = regsCombiGivenNumRegsFreq{q}(j,2:end) ;            
            
            idxNumRegCombi = idxNumSelRegsQ( ismember(selRegsHistoryGivenNumRegs,combi,'rows') == 1 ) ;
            thresholdsGivenModel{q}{j,2} = regProfileThresholdHistory(idxNumRegCombi,combi) ;
            degGivenModel{q}{j,2} = degHistory(idxNumRegCombi) ;
            if ratesFreqOn
                birthratesGivenModel{q}{j,2} = birthRatesHistory(idxNumRegCombi) ;
                m0GivenModel{q}{j,2} = m0History(idxNumRegCombi,:) ;
                switchesGivenModel{q}{j,2} = switchesHistory(idxNumRegCombi,:) ;
            end
            if ~isempty(idxNumRegCombi)
                thresholdsGivenModelKsDensQJ = zeros(120,length(combi)) ;
                for i = 1:length(combi)
                    lev = linspace(minRegProfiles(combi(i)),maxRegProfiles(combi(i)),120 ) ;
                    thresholdsGivenModelKsDensQJ(:,i) = ksdensity(regProfileThresholdHistory(idxNumRegCombi,combi(i)),lev) ;
                end
                thresholdsGivenModelKsDens{q}{j,2} = thresholdsGivenModelKsDensQJ ;
                lev = linspace(minDegRate,maxDegRate,120) ;
                degGivenModelKsDens{q}{j,2} = ksdensity(degHistory(idxNumRegCombi),lev) ;
                precGivenModel{q}{j,2} = precHistory(idxNumRegCombi,:) ;
                lsquaresWeightsPhiGivenModel{q}{j,2} = lsquaresWeightsPhiHistory(idxNumRegCombi,:) ;
                lev = linspace(0,12,120) ;
                for x = 1:numExps
                    precGivenModelKsDens{q}{j,2}(:,x) = ksdensity(precHistory(idxNumRegCombi,x),lev) ;
                end
                idxNumRegCombiNotEmpty(j)=1;
                thresholdsGivenModel{q}{j,1} =  combi;
                thresholdsGivenModelKsDens{q}{j,1} = combi ;
                degGivenModel{q}{j,1} = combi ;
                degGivenModelKsDens{q}{j,1} = combi ;
                precGivenModel{q}{j,1} = combi ;
                precGivenModelKsDens{q}{j,1} = combi ;
                lsquaresWeightsPhiGivenModel{q}{j,1} = combi ;
                birthratesGivenModel{q}{j,1} = combi ;
                m0GivenModel{q}{j,1} = combi ;
                switchesGivenModel{q}{j,1} = combi ;                
            end
            
            
        end
        thresholdsGivenModel{q} = thresholdsGivenModel{q}(logical(idxNumRegCombiNotEmpty),:);
        thresholdsGivenModelKsDens{q} = thresholdsGivenModelKsDens{q}(logical(idxNumRegCombiNotEmpty),:) ;
        degGivenModel{q} = degGivenModel{q}(logical(idxNumRegCombiNotEmpty),:) ;
        degGivenModelKsDens{q} = degGivenModelKsDens{q}(logical(idxNumRegCombiNotEmpty),:) ;
        precGivenModel{q} = precGivenModel{q}(logical(idxNumRegCombiNotEmpty),:) ;
        precGivenModelKsDens{q} = precGivenModelKsDens{q}(logical(idxNumRegCombiNotEmpty),:) ;
        lsquaresWeightsPhiGivenModel{q} = lsquaresWeightsPhiGivenModel{q}(logical(idxNumRegCombiNotEmpty),:) ;
    end
    
    output.regsCombiGivenNumRegsFreq = regsCombiGivenNumRegsFreq ;
    output.regsCombiFreq = regsCombiFreq ;
    output.modelFreq = modelFreq ;
    
    output.thresholdsGivenModel = thresholdsGivenModel ;
    output.thresholdsGivenModelKsDens = thresholdsGivenModelKsDens ;
    output.degGivenModel = degGivenModel ;
    output.degGivenModelKsDens = degGivenModelKsDens ;
    output.precGivenModel = precGivenModel ;
    output.precGivenModelKsDens = precGivenModelKsDens ;
    output.birthratesGivenModel = birthratesGivenModel ;
    output.m0GivenModel = m0GivenModel ;
    
    output.lsquaresWeightsPhiGivenModel = lsquaresWeightsPhiGivenModel ;
end

%% print model freqs
disp('Posterior distributions')
disp('--------------------')
disp('Number of regulators:')
numRegsFreq1 = [max(0,1-sum(numRegsFreq)) ; numRegsFreq] ;
for i=1:numRegs+1
    disp([num2str(i-1) ':prob=' num2str(numRegsFreq1(i))]) 
end
disp('   ----   ')
regsCombiFreq1 = regsCombiFreq{1} ;
disp('Regulators:')
for i=1:numRegs
    disp([regNames{i} ':prob=' num2str(regsCombiFreq1(i))]) 
end
disp('   ----   ')
disp('----------------------------------------')

bestModelprobs = zeros(numRegs,1) ; 
bestModelprobIdxs = zeros(numRegs,1) ;
for i=1:min(numRegs,maxNumRegs)    
    [bestModelprobs(i),bestModelprobIdxs(i)] = max(modelFreq{i}(:,1)) ;
end
[bestModelprob,bestModelprobIdx] = max(bestModelprobs) ;
bestModelregs = modelFreq{bestModelprobIdx}(bestModelprobIdxs(bestModelprobIdx),2:end) ;
numRegsBestModel = length(bestModelregs) ;
disp('Most likely Models')
disp('--------------------')
disp(['number of regulators:' num2str(numRegsBestModel)])
disp(['regulator set:' num2str(bestModelregs)])
disp(['probability:' num2str(bestModelprob)])
disp('    ----     ')
moreModels = 0 ;
modelFreq1 = modelFreq ;
bestModelprobIdx1 = bestModelprobIdx ;
bestModelprobIdxs1 = bestModelprobIdxs ;
while moreModels
    modelFreq1{bestModelprobIdx1}(bestModelprobIdxs1(bestModelprobIdx1),1) = 0 ;
    bestModelprobs1 = zeros(numRegs,1) ;
    bestModelprobIdxs1 = zeros(numRegs,1) ;
    for i=1:min(numRegs,maxNumRegs)  
        [bestModelprobs1(i),bestModelprobIdxs1(i)] = max(modelFreq1{i}(:,1)) ;
    end
    [bestModelprob1,bestModelprobIdx1] = max(bestModelprobs1) ;
    bestModelregs1 = modelFreq1{bestModelprobIdx1}(bestModelprobIdxs1(bestModelprobIdx1),2:end) ;
    numRegsBestModel1 = length(bestModelregs1) ;
    disp(['number of regulators:' num2str(numRegsBestModel1)])
    disp(['regulator set:' num2str(bestModelregs1)])
    disp(['probability:' num2str(bestModelprob1)])
    disp('    ----     ')
    prompt = 'display more models? (0/1)' ;
    if ~input(prompt)
        moreModels = 0 ;
    end
end
    

%% plots
if finalPlotsOn
    fs = 8 ;
    % get all regs to plot
    prompt = 'plot results (Y=1,N=0)?' ;
    if input(prompt)
        prompt = 'which regs to plot (vector)?' ;
        regsToPlot = input(prompt) ;
        numRegsToPlot = length(regsToPlot) ;
    else
        return
    end
    
    % plot all regs and posteriors
    for xxx=1
        targetmRNADataPointsColAll = cell(numExps,1) ;
        resAll = cell(numExps,1) ;
        
        q = numRegsToPlot ;
        combi = regsToPlot ;
        h = figure('Name','ProfileThreshActFit','NumberTitle','off','Position',[85.5839    1.1994   20.9903   25.6469]);
        set(h,'Units','centimeters')
        set(h,'Position',[2    4   (numExps+1)*7   (q+2)*3.2])
        
        lmargin = 0.07 ;
        bmarginw = 0.04 ;
        rmargin = 0.015 ;
        plotWidth = (1 - lmargin - numExps*bmarginw - rmargin)/(numExps+1) ;
        umargin = 0.05 ;
        bmarginh = 0.05 ;
        dmargin = 0.05 ;
        plotHeight = (1 - umargin - (q+1)*bmarginh - dmargin)/(q +2) ;
        for x = 1:numExps
            
            experiment = experiments{x,1};
            timepoints = experiment.timepoints; % observation times
            timepoints = timepoints - timepoints(1); % set to have first obs at time 0
            time = linspace(0,timepoints(end),timeFuncSize);
            
            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            plot(repmat(timepointsAll{x},nreplicatesAll(x),1),targetDataAll{x},'ok');
            if (x==1)
                ylabel('Target','FontSize',fs+2)
            end
            set(gca,'XTickMode','auto', 'FontSize',fs)
            axis([0 timepointsAll{x}(end) 0 maxTargetData+0.01*(maxTargetData - minTargetData)])
            hold on;
            if x>1
                set(gca,'YTickLabel',{})
            end
            
            figure(h);
            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            plot(time,ones(length(time),1),'w');
            hold on
            if (x==1)
                ylabel('\tau','FontSize',fs+3)
            end
            set(gca,'XTickMode','auto', 'FontSize',fs)
            axis('tight')
            if x>1
                set(gca,'YTickLabel',{})
            end
            set(gca,'XTickLabel',{})
            
            regProfileModelX = regProfileAll{x}(:,combi) ;
            thresholdsGivenModelKsDensQJ = zeros(profileFuncSize,length(combi)) ;
            for z = 1:q
                figure(h);
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                plot(time,regProfileModelX(:,z),'-k','linewidth',2) ;
                
                if (x == 1)
                    ylabel(regNames(combi(z)),'FontSize',fs+2)
                end
                if (z == 1)
                    title([experiment.name],'FontSize',fs+2)
                end
                set(gca,'XTickMode','auto', 'FontSize',fs)
                axis([0 timepoints(end) minRegProfiles(combi(z)) maxRegProfiles(combi(z))])
                if x>1
                    set(gca,'YTickLabel',{})
                end
                set(gca,'XTickLabel',{})
                
            end
            
        end
        for z = 1:q
            positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            hold on
            if (z == 1)
                title('Threshold posterior','FontSize',fs+2)
            end
            set(gca,'XTickMode','auto', 'FontSize',fs)
            set(gca,'YLim',[minRegProfiles(combi(z)) maxRegProfiles(combi(z))])
            set(gca,'YTickLabel',{})
            set(gca,'XTickLabel',{})
        end
        
        figure(h);
        positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight-dmargin/2] ;
        subplot('Position',positionVector)
        bar(0:numRegs,[1-sum(numRegsFreq) ; numRegsFreq],'k')
        set(gca,'XTick',0:numRegs)
        title('# of TFs','Fontsize',fs+2)
        set(gca,'XTickMode','auto', 'FontSize',fs)
        axis tight
        
        figure(h);
        positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin plotWidth  plotHeight-dmargin/2] ;
        subplot('Position',positionVector)
        bar(1:numRegs,regsCombiFreq{1}(:,1),'k')
        set(gca,'XTickMode','auto', 'FontSize',fs)
        title('which TF?','Fontsize',fs+2)
        set(gca,'XTick',1:numRegs)
        set(gca,'XTickLabel',regNames,'XTickLabelRotation',60)
        axis tight
    end
    
    % plot best model
    for xxx=1
        prompt = 'which regs to include in the model?' ;
        regsModel1 = input(prompt) ;
        numRegsModel1 = length(regsModel1) ;
        
        combi1 = regsModel1 ;
        q1 = numRegsModel1 ;
        
        if q1>q || sum(ismember(combi1,combi))~=length(combi1)
            prompt = 'select from the regulators included in the figure. new regulator set:';
            regsNew = input(prompt) ;
            numParModelNew = length(regsNew) ;
            q1 = numParModelNew ;
            combi1=regsNew;
        end
        
        combi1index = find(ismember(regsCombiGivenNumRegsFreq{q1}(:,2:end),combi1,'rows')) ;
        if interactivePlots
            while (regsCombiGivenNumRegsFreq{q1}(combi1index,1) == 0)
                disp('the selected regulator set is not sampled. Please select another set.');
                disp('Do you wish to see the posterior distributions for the sets of the same or smaller cardinality?')
                if input(prompt)
                    for i=1:q
                        disp(regsCombiGivenNumRegsFreq{i})
                    end
                end
                prompt = 'plot results? (0/1)' ;
                if input(prompt)
                    prompt = 'which regulators?' ;
                    regsNew = input(prompt) ;
                    numParModelNew = length(regsNew) ;
                    q1 = numParModelNew ;
                    combi1=regsNew;
                    combi1index = find(ismember(regsCombiGivenNumRegsFreq{q1}(:,2:end),combi1,'rows')) ;
                else
                    return
                end
            end
        end
        
        figure(h);
        positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight-dmargin/2] ;
        subplot('Position',positionVector)
        hold on
        y = zeros(numRegs+1,1) ;
        y(q1) = numRegsFreq(q1) ;
        bar(1:numRegs+1,y,'r')
        axis tight
        hold off
        
        figure(h);
        positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin plotWidth  plotHeight-dmargin/2] ;
        subplot('Position',positionVector)
        hold on
        y = zeros(numRegs,1) ;
        y(combi1) = regsCombiFreq{1}(combi1,1) ;
        bar(1:numRegs,y,'r')
        axis tight
        hold off
        
        % derive the best model for this combi
        for xxxx=1
            idxNumSelRegsQ = find(numSelRegsHistory == q1) ; % find which iterations had numRegs=q
            leNumSelRegsQ = length(idxNumSelRegsQ) ;  % how many they are
            selRegsHistoryGivenNumRegs = zeros(leNumSelRegsQ,q1) ;  % create a matrix with elements the regs for all iterations with numRegs=q
            for i = 1:leNumSelRegsQ
                selRegsHistoryGivenNumRegs(i,:) = selRegsHistory{idxNumSelRegsQ(i)}  ;
            end
            idxNumRegCombi = idxNumSelRegsQ( ismember(selRegsHistoryGivenNumRegs,combi1,'rows') == 1 ) ;
            
            levThresh = nan(profileFuncSize,q1) ;
%             thresholdsGivenModelQ = regProfileThresholdHistory(idxNumParParCombi,combi1) ;
            if ~isempty(regProfileThresholdHistory(idxNumRegCombi,combi1))
                thresholdsGivenModelKsDensQJ = zeros(profileFuncSize,q1) ; %thresholdsGivenModelKsDens{q}{j,2} ;
                for i = 1:q1
                    mini = min(regProfileThresholdHistory(idxNumRegCombi,combi1(i))) ;
                    maxi = max(regProfileThresholdHistory(idxNumRegCombi,combi1(i))) ;
%                     ri = maxi - mini ;
                    mini = max(mini,minRegProfiles(combi1(i))) ;
                    maxi = min(maxi,maxRegProfiles(combi1(i))) ;
                    levThresh(:,i) = linspace(mini,maxi,profileFuncSize) ;
                    thresholdsGivenModelKsDensQJ(:,i) = ksdensity(regProfileThresholdHistory(idxNumRegCombi,combi1(i)),levThresh(:,i)) ;%
                end
            end
            [~,modeidx] = max(thresholdsGivenModelKsDensQJ) ;
            th = zeros(q1,1) ;
            for i = 1:q1
                th(i) = levThresh(modeidx(i),i) ;
            end
%             th = median(thresholdsGivenModelQ) ;
            
            degGivenModelQ = degHistory(idxNumRegCombi) ;
            deg = median(degGivenModelQ) ;
            
            precGivenModelQ = precHistory(idxNumRegCombi,:) ;
            precGivenModelKsDensQ = zeros(120,numExps) ;
            for x=1:numExps
                lev = linspace(0,12,120) ;
                precGivenModelKsDensQ(:,x) = ksdensity(precGivenModelQ(:,x),lev) ;
            end
            precMeanGivenModelQ = mean(precGivenModelQ) ;
            
            lsquares = mcmcOutput.lsquares ;
            lsquares.weightsPhi = zeros(numExps,1) ;
            qAll = [] ;
            for x=1:numExps
                qx = repmat(lsquares.weights{x}.^(-lsquares.weightsPhi(x)),nreplicatesAll(x),1) ;
                qAll = [qAll ; qx]   ;
                lsquares.Q{x} = diag(qx) ;
            end
            lsquares.QAll = diag(qAll) ;
            
            switchesAll = cell(numExps,1) ;
            statesAll = cell(numExps,1) ;
            for x = 1:numExps
                
                experiment = experiments{x,1};
                timepoints = experiment.timepoints; % observation times
                timepoints = timepoints - timepoints(1); % set to have first obs at time 0
                time = linspace(0,timepoints(end),timeFuncSize);
                
                regProfileModelX = regProfileAll{x}(:,combi1) ;
                regProfileActivationX = zeros(timeFuncSize,q1);
                
                for z = 1:q1
                    %%% get activation
                    regProfileActivationX(hSteps(x):end,z) = getRegProfileActivation(regProfileModelX(hSteps(x):end,z),th(z));
                    % set the act func before delay as the first after delay to avoid switches
                    regProfileActivationX(1:hSteps(x)-1,z) = regProfileActivationX(hSteps(x),z) ;
                end
                %%% get vector with state at each (dense) time point for all regs
                twosMat = repmat(pow2(0:(q1-1)),timeFuncSize,1);
                stateChanges = getGeneActivation(regProfileActivationX,twosMat);
                %%% get states used and switch times
                %%% if too fast switches, get them to the largest switch and set
                %%% the state as after the last
                [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x)); % states used and at which times
                switchesAll{x} = switches;
                statesAll{x} = states;
            end
            %%% get the states used across experiments
            statesUsed = getStatesUsed(statesAll);
            
            %%%% perfrom multi-exp regression
            %%% get m0's and tau-rates
            YAll = profilesOutput.YAll ;
            [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,nreplicatesAll,lsquares) ;
            
            %%% for wls I'll get the above (ols) results, re-compute phi and
            %%% perform wls
            if strcmp(lsquares.type,'wls')
                e_hat = cell(numExps,1) ;
                targetmRNA = cell(numExps,1) ;
                qAll = [] ;
                for x = 1:numExps
                    
                    %         keyboard
                    br = zeros(length(statesAll{x}),1) ;
                    for y=1:length(br)
                        br(y) = find(statesAll{x}(y) == statesUsed) ;
                    end
                    targetmRNA{x} = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                    
                    %%% calculate SSE of profile to data
                    e_hat{x} = calculateEHat(targetmRNA{x},targetDataAll{x},eHatMatAll{x},nreplicatesAll(x)) ;
                    
                    yphi = log(precMeanGivenModelQ(x)*e_hat{x}.^2) ;
                    xphi = -repmat(log(lsquares.weights{x}),nreplicatesAll(x),1) ;
                    phi_hat = (xphi'*xphi)^(-1)*xphi'*yphi ;
                    phi_hat = max(min(phi_hat,1),0);
                    lsquares.weightsPhi(x) = phi_hat ;
                    qx = repmat(lsquares.weights{x}.^(-lsquares.weightsPhi(x)),nreplicatesAll(x),1) ;
                    qAll = [qAll ; qx]   ;
                    lsquares.Q{x} = diag(qx) ;
                end
                QAll = diag(qAll) ;
                lsquares.QAll = QAll ;
                
                [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,nreplicatesAll,lsquares) ;
                
            end
        end
        
        rmin = [] ;
        rmax = [] ;
        for x = 1:numExps
            
            experiment = experiments{x,1};
            timepoints = experiment.timepoints; % observation times
            timepoints = timepoints - timepoints(1); % set to have first obs at time 0
            time = linspace(0,timepoints(end),timeFuncSize);
            
            br = zeros(length(statesAll{x}),1) ;
            for y=1:length(br)
                br(y) = find(statesAll{x}(y) == statesUsed) ;
            end
            targetmRNA = nStateSwitchODE(time,m0All(x),birthRates(br),deg,switchesAll{x});
            
            targetmRNADataPoints = nStateSwitchODE(timepoints,m0All(x),birthRates(br),deg,switchesAll{x});
            
            targetmRNADataPointsCol = repmat(targetmRNADataPoints,1,experiment.nreplicates)' ;
            childDataAllCol = reshape(targetDataAll{x},experiment.nreplicates*experiment.ntimepoints,1);
            res = childDataAllCol - targetmRNADataPointsCol ;
            
            targetmRNADataPointsColAll{x} = targetmRNADataPointsCol ;
            resAll{x} = res ;
            
            t = 0;
            r = birthRates(br(1));
            for z = 1:length(switchesAll{x})
                t = [t switchesAll{x}(z) switchesAll{x}(z)];
                r = [r birthRates(br(z)) birthRates(br(z+1))];
            end
            t = [t time(end)];
            r = [r birthRates(br(end))];
            
            rmin = min([r rmin]) ;
            rmax = max([r rmax]) ;
            
            figure(h);
            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            hold on
            plot(t,r,'-r','linewidth',2);
            axis([t(1) t(end) rmin-0.05*(rmax-rmin)-1 rmax+0.05*(rmax-rmin)+1])
            
            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            hold on
            plot(time,targetmRNA,'-r','linewidth',2);
            
            for z1 = 1:q1
                z = find(combi==combi1(z1));
                figure(h);
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                hold on
                plot([timepoints(1) timepoints(end)+timepoints(2)],th(z1)*ones(2,1),'r--','linewidth',2)
                set(gca,'XTickMode','auto', 'FontSize',fs)
                axis([0 timepoints(end) minRegProfiles(combi1(z1)) maxRegProfiles(combi1(z1))])                                
            end
            
        end
        for z1 = 1:q1
            z = find(combi==combi1(z1));
            figure(h);
            positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            hold on
            plot(thresholdsGivenModelKsDensQJ(:,z1),levThresh(:,z1),'-xr','linewidth',2,'markersize',0.25)
            set(gca,'YLim',[minRegProfiles(combi1(z1)) maxRegProfiles(combi1(z1))])
            set(gca,'YTickLabel',{})
            xlim = get(gca,'XLim') ;
            if z1==q
                set(gca,'XTickLabel',{num2str(round(xlim(1),3)); num2str(round(xlim(end)/2,3)); num2str(round(xlim(end),3))})
            else
                set(gca,'XTickLabel',{})
            end
        end
    end
    
    % more models
    if plot2modelsOn
        for xxxxx=1
            % get model and check if sampled
            for xxxx=1
                
                prompt = 'plotting other models? 1=Y/0=N:' ;
                if input(prompt)
                    prompt = 'which set of regs (vector)?' ;
                    regsNew = input(prompt) ;
                    numParModelNew = length(regsNew) ;
                    q2 = numParModelNew ;
                    combi2=regsNew;
                    if ~isempty(setdiff(regsNew,regsToPlot))
                        disp('some of the selected parents are not included in the current figure')
                        return
                    end
                else
                    return
                end
                
                
                if q2>q || sum(ismember(combi1,combi))~=length(combi1)
                    prompt = 'select from the regulators included in the figure. new regulator set:';
                    regsNew = input(prompt) ;
                    numParModelNew = length(regsNew) ;
                    q2 = numParModelNew ;
                    combi2=regsNew;
                end
                
                % check if the selected model is sampled
                combi2index = find(ismember(regsCombiGivenNumRegsFreq{q2}(:,2:end),combi2,'rows')) ;
                if interactivePlots
                    while (regsCombiGivenNumRegsFreq{q2}(combi2index,1) == 0)
                        disp('the selected regulator set is not sampled. Please select another set.');
                        disp('The posterior distributions for the sets of the same or smaller cardinality are')
                        for i=1:q
                            disp(regsCombiGivenNumRegsFreq{i})
                        end
                        prompt = 'plot results? (0/1)' ;
                        if input(prompt)
                            prompt = 'which regulators?' ;
                            regsNew = input(prompt) ;
                            numParModelNew = length(regsNew) ;
                            q2 = numParModelNew ;
                            combi2=regsNew;
                            combi2index = find(ismember(regsCombiGivenNumRegsFreq{q2}(:,2:end),combi2,'rows')) ;
                        else
                            return
                        end
                    end
                end
            end
            
            % get model parameters
            for xxxx=1
                idxNumSelRegsQ = find(numSelRegsHistory == q2) ; % find which iterations had numRegs=q
                leNumSelRegsQ = length(idxNumSelRegsQ) ;  % how many they are
                selRegsHistoryGivenNumRegs = zeros(leNumSelRegsQ,q2) ;  % create a matrix with elements the regs for all iterations with numRegs=q
                for i = 1:leNumSelRegsQ
                    selRegsHistoryGivenNumRegs(i,:) = selRegsHistory{idxNumSelRegsQ(i)}  ;
                end
                
                idxNumRegCombi = idxNumSelRegsQ( ismember(selRegsHistoryGivenNumRegs,combi2,'rows') == 1 ) ;
%                 thresholdsGivenModelQ = regProfileThresholdHistory(idxNumParParCombi,combi2) ;
                if ~isempty(regProfileThresholdHistory(idxNumRegCombi,combi2))
                    thresholdsGivenModelKsDensQJ = zeros(profileFuncSize,q2) ;
                    levThresh = nan(profileFuncSize,q2) ;
                    for i = 1:q2
                        mini = min(regProfileThresholdHistory(idxNumRegCombi,combi2(i))) ;
                        maxi = max(regProfileThresholdHistory(idxNumRegCombi,combi2(i))) ;
%                         ri = maxi - mini ;
                        mini = max(mini,minRegProfiles(combi2(i))) ;
                        maxi = min(maxi,maxRegProfiles(combi2(i))) ;
                        levThresh(:,i) = linspace(mini,maxi,profileFuncSize) ;
                        thresholdsGivenModelKsDensQJ(:,i) = ksdensity(regProfileThresholdHistory(idxNumRegCombi,combi2(i)),levThresh(:,i)) ;
                    end
                end
                [~,modeidx] = max(thresholdsGivenModelKsDensQJ) ;
                th = zeros(q2,1) ;
                for i = 1:q2
                    th(i) = levThresh(modeidx(i),i) ;
                end
%                 th = median(thresholdsGivenModelQ) ;
                
                degGivenModelQ = degHistory(idxNumRegCombi) ;
                deg = median(degGivenModelQ) ;
                
                precGivenModelQ = precHistory(idxNumRegCombi,:) ;
                precGivenModelKsDensQ = zeros(120,numExps) ;
                for x=1:numExps
                    lev = linspace(0,12,120) ;
                    precGivenModelKsDensQ(:,x) = ksdensity(precGivenModelQ(:,x),lev) ;
                end
                precMeanGivenModelQ = mean(precGivenModelQ);
                
                lsquares = mcmcOutput.lsquares ;
                lsquares.weightsPhi = zeros(numExps,1) ;
                qAll = [] ;
                for x=1:numExps
                    qx = repmat(lsquares.weights{x}.^(-lsquares.weightsPhi(x)),nreplicatesAll(x),1) ;
                    qAll = [qAll ; qx]   ;
                    lsquares.Q{x} = diag(qx) ;
                end
                lsquares.QAll = diag(qAll) ;
                
                switchesAll = cell(numExps,1) ;
                statesAll = cell(numExps,1) ;
                for x = 1:numExps
                    
                    experiment = experiments{x,1};
                    timepoints = experiment.timepoints; % observation times
                    timepoints = timepoints - timepoints(1); % set to have first obs at time 0
                    time = linspace(0,timepoints(end),timeFuncSize);
                    
                    regProfileModelX = regProfileAll{x}(:,combi2) ;
                    regProfileActivationX = zeros(timeFuncSize,q2);
                    
                    for z = 1:q2
                        %%% get activation
                        regProfileActivationX(hSteps(x):end,z) = getRegProfileActivation(regProfileModelX(hSteps(x):end,z),th(z));
                        % set the act func before delay as the first after delay to avoid switches
                        regProfileActivationX(1:hSteps(x)-1,z) = regProfileActivationX(hSteps(x),z) ;
                    end
                    
                    %%% get vector with state at each (dense) time point for all regs
                    twosMat = repmat(pow2(0:(q2-1)),timeFuncSize,1);
                    stateChanges = getGeneActivation(regProfileActivationX,twosMat);
                    %%% get states used and switch times
                    %%% if too fast switches, get them to the largest switch and set
                    %%% the state as after the last
                    [switches,states] = getSwitchPoints(stateChanges,time,minSwTime(x)); % states used and at which times
                    switchesAll{x} = switches;
                    statesAll{x} = states;
                end
                %%% get the states used across experiments
                statesUsed = getStatesUsed(statesAll);
                
                %%%% perfrom multi-exp regression
                %%% get m0's and tau-rates
                YAll = profilesOutput.YAll ;
                [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,nreplicatesAll,lsquares) ;
                
                %%% for wls I'll get the above (ols) results, re-compute phi and
                %%% perform wls
                if strcmp(lsquares.type,'wls')
                    e_hat = cell(numExps,1) ;
                    targetmRNA = cell(numExps,1) ;
                    qAll = [] ;
                    for x = 1:numExps
                        br = zeros(length(statesAll{x}),1) ;
                        for y=1:length(br)
                            br(y) = find(statesAll{x}(y) == statesUsed) ;
                        end
                        targetmRNA{x} = nStateSwitchODE(timepointsAll{x},m0All(x),birthRates(br),deg,switchesAll{x});
                        
                        %%% calculate SSE of profile to data
                        e_hat{x} = calculateEHat(targetmRNA{x},targetDataAll{x},eHatMatAll{x},nreplicatesAll(x)) ;
                        
                        yphi = log(precMeanGivenModelQ(x)*e_hat{x}.^2) ;
                        xphi = -repmat(log(lsquares.weights{x}),nreplicatesAll(x),1) ;
                        phi_hat = (xphi'*xphi)^(-1)*xphi'*yphi ;
                        phi_hat = max(min(phi_hat,1),0);
                        lsquares.weightsPhi(x) = phi_hat ;
                        qx = repmat(lsquares.weights{x}.^(-lsquares.weightsPhi(x)),nreplicatesAll(x),1) ;
                        qAll = [qAll ; qx]   ;
                        lsquares.Q{x} = diag(qx) ;
                    end
                    QAll = diag(qAll) ;
                    lsquares.QAll = QAll ;
                    
                    [m0All,birthRates] = calculateFixedStateModel(timepointsAll,switchesAll,statesAll,statesUsed,deg,YAll,ntimepointsAll,nreplicatesAll,lsquares) ;
                    
                end
                
            end
            
            for x = 1:numExps
                
                experiment = experiments{x,1};
                timepoints = experiment.timepoints; % observation times
                timepoints = timepoints - timepoints(1); % set to have first obs at time 0
                time = linspace(0,timepoints(end),timeFuncSize);
                
                br = zeros(length(statesAll{x}),1) ;
                for y=1:length(br)
                    br(y) = find(statesAll{x}(y) == statesUsed) ;
                end
                targetmRNA = nStateSwitchODE(time,m0All(x),birthRates(br),deg,switchesAll{x});
                
                targetmRNADataPoints = nStateSwitchODE(timepoints,m0All(x),birthRates(br),deg,switchesAll{x});
                
                targetmRNADataPointsCol = repmat(targetmRNADataPoints,1,experiment.nreplicates)' ;
                childDataAllCol = reshape(targetDataAll{x},experiment.nreplicates*experiment.ntimepoints,1);
                res = childDataAllCol - targetmRNADataPointsCol ;
                
                targetmRNADataPointsColAll{x} = targetmRNADataPointsCol ;
                resAll{x} = res ;
                
                t = 0;
                r = birthRates(br(1));
                for z = 1:length(switchesAll{x})
                    t = [t switchesAll{x}(z) switchesAll{x}(z)];
                    r = [r birthRates(br(z)) birthRates(br(z+1))];
                end
                t = [t time(end)];
                r = [r birthRates(br(end))];
                
                rmin = min([r rmin]) ;
                rmax = max([r rmax]) ;
                
                figure(h);
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight] ;
                subplot('Position',positionVector) ;
                hold on
                plot(t,r,'-c','linewidth',2);
                axis([t(1) t(end) rmin-0.05*(rmax-rmin)-1 rmax+0.05*(rmax-rmin)+1])
                
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                hold on;
                plot(time,targetmRNA,'-c','linewidth',2);
                
                for z2 = 1:q2                    
                    z = find(combi == combi2(z2)) ;
                    figure(h);
                    positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
                    subplot('Position',positionVector)
                    hold on
                    plot(time,th(z2)*ones(length(time),1),'--c','linewidth',2)
                    set(gca,'XTickMode','auto', 'FontSize',fs)
                    axis([0 timepoints(end) minRegProfiles(combi2(z2)) maxRegProfiles(combi2(z2))])                                        
                end
                
            end
            for z2 = 1:q2
                z = find(combi == combi2(z2)) ;
                positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+(2+q-z)*(bmarginh+plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)%
                hold on                
                plot(thresholdsGivenModelKsDensQJ(:,z2),levThresh(:,z2),'-c','linewidth',2)
                set(gca,'YLim',[minRegProfiles(combi2(z2)) maxRegProfiles(combi2(z2))])
                set(gca,'YTickLabel',{})
                set(gca,'XTickLabel',{})
            end
            
            figure(h);
            positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin+plotHeight+bmarginh plotWidth  plotHeight-dmargin/2] ;
            subplot('Position',positionVector)
            hold on
            bar(q2,numRegsFreq(q2),'c')
            
            figure(h);
            positionVector = [lmargin + numExps*(bmarginw + plotWidth) dmargin plotWidth  plotHeight-dmargin/2] ;
            subplot('Position',positionVector)
            hold on
            y = zeros(numRegs,1) ;
            y(combi2) = regsCombiFreq{1}(combi2,1) ;
            bar(1:numRegs,y,'c')
            hold off
            
        end
    end
end

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

function statesUsed = getStatesUsed(states)

statesUsed = states{1};

for x = 2:length(states)
    statesUsed = unique([states{x} statesUsed]);
end

end

function [e_hat] = calculateEHat(profile,data,eHatMat,nReps)
for x = 1:nReps
    eHatMat(x,:) = profile ; % create the fit matrix with nrows = nreplicates, ncolumns = timepoints
end
e_hat = eHatMat - data ;% compute the residuals
e_hatT = e_hat' ;
e_hat = e_hatT(:) ;
end


