function output = createPriorsProp(profilesOutput,priors,nregsPriors,regsPriors,threshPriors,...
    maxNumRegs,hyperparamsLambda0,...
    regulatorsPriorType,thresholdsPriorType,thresholdsPriorLBoundStandGrads,...
    hyperparamsDeg,hyperparamsMinDegRate,hyperparamsMaxDegRate,hyperparamsPrior_deg_nu0,...
    hyperparamsPrior_deg_s02,hyperparamsPrecDeg0,hyperparamsPrecSigma0,...
    profileFuncSize,priorPlotsOn,numRegulatorsInFig)

addpath(genpath(pwd));

%% get variables from input structures 
experiments = profilesOutput.experiments ;
regulatorNames = profilesOutput.regulatorNames ;

nexps = size(experiments,1);
nregulators = length(regulatorNames);

timeFuncSize = profilesOutput.timeFuncSize ;

%% collect various hyperparameters
%% deg
deg = hyperparamsDeg ;  
minDegRate = hyperparamsMinDegRate ; 
maxDegRate = hyperparamsMaxDegRate ; 
prior_deg_nu0 = hyperparamsPrior_deg_nu0 ; 
prior_deg_s02 = hyperparamsPrior_deg_s02 ; 
prior_deg_a0 = prior_deg_nu0/2 ;  % a,b hyperparameters for gamma prior of deg (equal to non-central chi2 distribution)
prior_deg_b0 = 2/(prior_deg_nu0*prior_deg_s02) ;
%% precision prior hyperparameters
precDeg0 = hyperparamsPrecDeg0 ;
precSigma0 = hyperparamsPrecSigma0 ;
%% protein ODE
trans_p = ones(nregulators,1);  % translation rate
deg_p = ones(nregulators,1) .* 0.075; % protein degradation
%% the max number of selected parents for a single target
%% mean (lambda) of the poisson prior for the # of selected parents
numSelRegsPriorLambda = hyperparamsLambda0 ;

%% compute nregulators prior
if nregsPriors
    nregulatorsPriorOutput = priors.nregulatorsPriorOutput ;
else
    nregulatorsPrior = zeros(nregulators+1,1) ;
    for z = 0:maxNumRegs
        nregulatorsPrior(z+1) = truncPoisson(z,numSelRegsPriorLambda,maxNumRegs) ;
    end
    nregulatorsPriorOutput.probs = nregulatorsPrior ;
    nregulatorsPriorOutput.type = 'poisson' ;
end

%% get parents prior 
if regsPriors
    regsPriorOutput = priors.regsPriorOutput ;
else
    if strcmp(regulatorsPriorType,'range')
        nSplineSamples = profilesOutput.nSplineSamples ;
        
        minRegProfileBoot = profilesOutput.minRegProfileBoot ;
        maxRegProfileBoot = profilesOutput.maxRegProfileBoot ;
        
        % compute range in each Bootstrap sample
        rangeReg = zeros(nSplineSamples,nregulators) ;
        for z=1:nregulators
            for n=1:nSplineSamples
                rangeReg(n,z) = maxRegProfileBoot(n,z) - minRegProfileBoot(n,z) ;
            end
        end
        
        regsPriorOutput = createRegulatorsPrior(rangeReg)   ;
        regsPriorOutput.type = 'range' ;
    else
        regsPriorOutput.type = 'uniform' ;
        regsPriorOutput.regPriorsMatrixAllcomb = [] ;
    end
end


%% get thresholds prior
if threshPriors
    thresholdsPriorOutput = priors.thresholdsPriorOutput ; 
else
    if strcmp(thresholdsPriorType,'gradient')
        lBoundStandGrads = thresholdsPriorLBoundStandGrads ;
        thresholdsPriorOutput = createThresholdsPrior(profilesOutput,timeFuncSize,profileFuncSize,lBoundStandGrads) ;
        thresholdsPriorOutput.type = 'gradient';
    else
        thresholdsPriorOutput.type = 'uniform' ;
        thresholdsPriorOutput.threshPriorDensity2 = [] ;
    end
end

%% get plots
if (strcmp(thresholdsPriorType,'gradient') && priorPlotsOn)
    counter = nregulators ;
    k = 0 ;
    while (counter > 0)
        m1 = min(nregulators-k,numRegulatorsInFig) ;
        figure('name','dataANDpriors','numbertitle','off','windowstyle','docked');
        lmargin = 0.1 ;
        bmarginw = 0.04 ;
        rmargin = 0.03 ;
        plotWidth = (1 - lmargin - nexps*bmarginw - rmargin)/(nexps+1) ;
        
        umargin = 0.05 ;
        bmarginh = 0.05 ;
        dmargin = 0.05 ;
        plotHeight = (1 - umargin - m1*bmarginh - dmargin)/(m1 +1) ;
        
        fs=8;
        timepointsAll = profilesOutput.timepointsAll ;
        nreplicatesAll = profilesOutput.nreplicatesAll ;
        targetDataAll = profilesOutput.targetDataAll ;
        timeAll = profilesOutput.timeAll ;
        regProfileAll = profilesOutput.regProfileAll ;
        minTargetData = profilesOutput.minTargetData ;
        maxTargetData = profilesOutput.maxTargetData ;
        regProfileAllBoot = profilesOutput.regProfileAllBoot   ;
        nregulators = profilesOutput.nregulators ;
        minRegProfiles = profilesOutput.minRegProfiles ;
        maxRegProfiles = profilesOutput.maxRegProfiles ;
        for x=1:nexps
            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            plot(repmat(timepointsAll{x},nreplicatesAll(x),1),targetDataAll{x},'or','linewidth',2);
            if (x == 1)
                ylabel('Target','FontSize',fs+2,'FontWeight','bold')
            end
            xlabel('time','Fontsize',fs)
            axis([0 timepointsAll{x}(end) 0 maxTargetData+0.1*(maxTargetData - minTargetData)])
            set(gca,'XTickMode','auto', 'FontSize',fs,'FontWeight','bold')
            xtick = get(gca,'XTick') ;
            if x>1
                set(gca,'YTick',[])
            end
            
            for z=k+1:k+m1
                
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin + (m1-(z-k)+1)*(bmarginh + plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                plot(timeAll{x},regProfileAllBoot{x}(:,1:10:end,z),'c-','linewidth',2)
                hold all
                plot(timeAll{x},regProfileAll{x}(:,z),'-k','linewidth',2)
                set(gca,'XTick',xtick)
                set(gca,'XtickLabel',[])
                if (x == 1)
                    ylabel(['reg' num2str(z)],'FontSize',fs,'FontWeight','bold')
                end
                if x>1
                    set(gca,'YtickLabel',[])
                end
                set(gca,'XTickMode','auto', 'FontSize',fs,'FontWeight','bold')
                if (z == 1)
                    title(['Experiment ' num2str(x)],'FontSize',fs,'FontWeight','bold')
                end
                axis([0 timepointsAll{x}(end) minRegProfiles(z)-0.01*(maxRegProfiles(z) - minRegProfiles(z)) maxRegProfiles(z)+0.01*(maxRegProfiles(z) - minRegProfiles(z))])
                %             axis([0 timescale(end) minRegProfilesBootAll(z)-0.05*(maxRegProfilesBootAll(z) - minRegProfilesBootAll(z)) maxRegProfilesBootAll(z)+0.05*(maxRegProfilesBootAll(z) - minRegProfilesBootAll(z))])
                yticks = get(gca,'YTick') ;
                if length(yticks)>2
                    set(gca,'YTick',yticks(1:2:end))
                end
                
                
            end
            
        end
        for z=k+1:k+m1
            le = linspace(minRegProfiles(z),maxRegProfiles(z),profileFuncSize) ;
            posX = lmargin + nexps*(bmarginw + plotWidth) ;
            posY = dmargin + (m1-(z-k)+1)*(bmarginh + plotHeight) ;
            positionVector = [posX posY plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            %         subplot(numRegents+1,numExps+1,(z-1)*(numExps+1) + (numExps+1))
            plot(thresholdsPriorOutput.threshPriorDensity2(:,z),le,'b-','linewidth',2)
            hold all
            axis normal
            set(gca,'YLim',[minRegProfiles(z)-0.01*(maxRegProfiles(z) - minRegProfiles(z)) maxRegProfiles(z)+0.01*(maxRegProfiles(z) - minRegProfiles(z))])
            if (z == 1)
                title('Threshold prior','FontSize',fs+2,'FontWeight','bold')
            end
            set(gca,'YTickLabel',{})
            set(gca,'XTickMode','auto', 'FontSize',fs,'FontWeight','bold')
            set(gca,'XtickLabel',{})
        end     
        counter = counter - m1 ;
        k = k + numRegulatorsInFig ;
    end
    positionVector = [lmargin + (nexps)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
    subplot('Position',positionVector)
    pp = regsPriorOutput.regPriorsMatrixAllcomb{1};
    bar(1:nregulators,pp,'b')
    xlabel('Regulator prior','Fontsize',fs)
    set(gca,'XTickMode','auto', 'FontSize',fs,'FontWeight','bold')
    set(gca,'XTick',1:nregulators)
    axis tight
end

% prepare output
for xx=1
    priors = struct ;
    priors.deg = deg ;
    priors.minDegRate = minDegRate ; % minimum value for the bounded proposal ration for deg_m
    priors.maxDegRate = maxDegRate ; % maximum value for the bounded proposal ration for deg_m
    priors.prior_deg_a0 = prior_deg_a0 ;
    priors.prior_deg_b0 = prior_deg_b0 ;
    
    priors.precDeg0 = precDeg0 ;
    priors.precSigma0 = precSigma0 ;
    
    priors.trans_p = trans_p ;
    priors.deg_p = deg_p ;
    priors.numSelRegsPriorLambda = numSelRegsPriorLambda ;
    
    priors.nregulatorsPriorOutput = nregulatorsPriorOutput ;
    priors.regsPriorOutput = regsPriorOutput ;
    priors.thresholdsPriorOutput = thresholdsPriorOutput ;
    
    output = struct;
    output.priors = priors ;
    output.profileFuncSize = profileFuncSize ;
    output.maxNumRegs = maxNumRegs ;
end

end


function prob = truncPoisson(k,lambda,kmax)

if (k > -1 && k < kmax+1)
    prob = exp(-lambda)*lambda^k/factorial(k) ;
else
    prob = 0 ;
end

end
