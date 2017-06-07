function output = getProfiles(targetName,regulatorNames,datasets,smoothProfiles,regProfilePlots,numRegulatorsInFig,minSwTimeDefault,minSwTime)

% gets the profiles of regulators and target data
% for regulators data are collected from the createSmoothProfiles output


nexps = size(datasets,1);
nregulators = size(regulatorNames,1);

%% size of functions created via smoothing etc.
timeFuncSize = 100;
nBootSamples = smoothProfiles.nBootSamples ;
hSteps = smoothProfiles.hSteps ;

%% create all variables
timepointsAll = cell(nexps,1);
timeAll = cell(nexps,1);
timepointsColAll = cell(nexps,1) ;
expNameAll = cell(nexps,1);
ntimepointsAll = zeros(nexps,1);
nreplicatesAll = zeros(nexps,1);

%% regulators expression, activation, etc.
regProfileAll = cell(nexps,1) ;
regProfileGradAll = cell(nexps,1) ;
regProfileAbsGradAll = cell(nexps,1) ;
regProfileAllBoot = cell(nexps,1);
regProfileGradAllBoot = cell(nexps,1);
regProfileAbsGradAllBoot = cell(nexps,1);
regProfileActivationAll = cell(nexps,1);
regProfileActivationAllProposed = cell(nexps,1);
eHatMatAll = cell(nexps,1);
sampleSizeAll = zeros(nexps,1) ;
targetIDAll = zeros(nexps,1);
targetDataAll = cell(nexps,1);
targetSmoothData = cell(nexps,1) ;
if minSwTimeDefault
    minSwTime = zeros(nexps,1);
end

%% various objects used in the algorithm
currIdx = 1; % used to construct the responses vector
minTargetData = [] ;
maxTargetData = [] ;
for x = 1:nexps
    experiment = datasets{x,1};
    
    targetID = findIndex(targetName,experiment.orderedNames);
    targetIDAll(x) = targetID;
    
    expName = experiment.name; 
    ntimepoints = experiment.ntimepoints;
    nreplicates = experiment.nreplicates; 
    
    expNameAll{x} = expName;
    ntimepointsAll(x) = ntimepoints;
    nreplicatesAll(x) = nreplicates;
        
    sampleSizeAll(x) = ntimepoints*nreplicates ;
    
    d = reshape(experiment.data(targetID,:),ntimepoints,nreplicates);
    YAll(currIdx:(currIdx - 1 + ntimepoints * nreplicates)) = reshape(d,ntimepoints * nreplicates,1);
    currIdx = currIdx + (ntimepoints * nreplicates);
    data = d';
    targetDataAll{x} = data;
    
    minTargetData = min([minTargetData min(data)]) ;
    maxTargetData = max([maxTargetData max(data)]) ;
    
    timepoints = experiment.timepoints; 
    timepoints = timepoints - timepoints(1); 
    timepointsAll{x} = timepoints;
    
    timepointsCol = repmat(timepoints,1,nreplicates)' ;
    timepointsColAll{x} = timepointsCol ;
    
    targetDataCol = reshape(data',1,ntimepoints*nreplicates) ;
    
    ftarget = fit(timepointsCol,targetDataCol','smoothingspline');
    targetSmoothData{x} = feval(ftarget,timepoints);
    
    %%% min time between switches
    if minSwTimeDefault
        deltaTime = min(timepoints(2:end)-timepoints(1:(end-1)));
        minSwTime(x) = deltaTime; 
    end
    
    time = linspace(0,timepoints(end),timeFuncSize);
    timeAll{x} = time;
    
end
YAll = YAll';

regNamesAll = smoothProfiles.regulatorNames ;
regidxs = nan(nregulators,1) ;
for z=1:nregulators
    regidxs(z) = find(strcmp(regulatorNames(z),regNamesAll)) ;
end
        
regDataAll = smoothProfiles.regDataAll(:,regidxs) ;
regDataColAll = smoothProfiles.regDataColAll(:,regidxs) ;
smoothRegProfileAllBoot = smoothProfiles.smoothRegProfileAllBoot(:,regidxs);
smoothRegProfileAll = smoothProfiles.smoothRegProfileAll(:,regidxs);
regIDAll = smoothProfiles.regIDAll(:,regidxs);

for x=1:nexps
    regProfileAll{x} = smoothProfiles.regProfileAll{x}(:,regidxs) ;
    regProfileGradAll{x} = smoothProfiles.regProfileGradAll{x}(:,regidxs) ;
    regProfileAbsGradAll{x} = smoothProfiles.regProfileAbsGradAll{x}(:,regidxs) ;
    regProfileAllBoot{x} = smoothProfiles.regProfileAllBoot{x}(:,:,regidxs);
    regProfileGradAllBoot{x} = smoothProfiles.regProfileGradAllBoot{x}(:,:,regidxs);
    regProfileAbsGradAllBoot{x} = smoothProfiles.regProfileAbsGradAllBoot{x}(:,:,regidxs);
    
    regProfileActivationAll{x} = smoothProfiles.regProfileActivationAll{x}(:,regidxs);
    regProfileActivationAllProposed{x} =  smoothProfiles.regProfileActivationAllProposed{x}(:,regidxs);
end

minRegProfileBoot = smoothProfiles.minRegProfilesBoot(:,regidxs) ;
maxRegProfileBoot = smoothProfiles.maxRegProfilesBoot(:,regidxs) ;

minRegProfiles = smoothProfiles.minRegProfiles(regidxs) ;
maxRegProfiles = smoothProfiles.maxRegProfiles(regidxs) ;

minDataReg = smoothProfiles.minDataReg(regidxs) ;
maxDataReg = smoothProfiles.maxDataReg(regidxs) ;

minData = min(minDataReg) ;
maxData = max(maxDataReg) ;

resAll = smoothProfiles.resAll(:,regidxs) ;
smoothRegProfileColAll = smoothProfiles.smoothRegProfileColAll(:,regidxs) ;


counter = nregulators ;
k = 0 ;
if regProfilePlots
    while (counter > 0)
        m1 = min(nregulators-k,numRegulatorsInFig) ;
        figure('name','data','numbertitle','off','windowstyle','docked');
        lmargin = 0.1 ;
        bmarginw = 0.04 ;
        rmargin = 0.03 ;
        plotWidth = (1 - lmargin - nexps*bmarginw - rmargin)/nexps ;
       
        umargin = 0.05 ;
        bmarginh = 0.05 ;
        dmargin = 0.05 ;
        plotHeight = (1 - umargin - m1*bmarginh - dmargin)/(m1 +1) ;
        
        for x=1:nexps
            experiment = datasets{x};
            expName = expNameAll{x};
            timepoints = experiment.timepoints; % observation times
            timepoints = timepoints - timepoints(1); % set to have first obs at time 0
            nreplicates = experiment.nreplicates; % number of replicates

            time = linspace(0,timepoints(end),timeFuncSize);

            positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin  plotWidth  plotHeight] ;
            subplot('Position',positionVector)
            plot(repmat(timepoints,nreplicates,1),targetDataAll{x},'or');
            if (x == 1)
                ylabel('target','FontSize',12,'FontWeight','bold')
            end
            axis([0 timepoints(end) minTargetData-0.01*(maxTargetData - minTargetData) maxTargetData+0.01*(maxTargetData - minTargetData)])
            
            for z=k+1:k+m1
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin + (m1-(z-k)+1)*(bmarginh + plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                plot(repmat(timepoints,1,nreplicates),regDataColAll{x,z},'ob');
                hold all
                plot(time,smoothRegProfileAllBoot{x,z}(:,1),'-k')
                if (z == 1)
                    title(['Experiment ' expName],'FontSize',12.5,'FontWeight','bold')
                end
                if (x == 1)
                    ylabel(['reg ' num2str(z)],'FontSize',12,'FontWeight','bold')
                end
                
                for n=2:10:nBootSamples
                    plot(time,smoothRegProfileAllBoot{x,z}(:,n),'-c')
                    hold all
                end
                plot(time,smoothRegProfileAll{x,z},'-.k','lineWidth',1.5)
                axis tight
            end
        end
        counter = counter - m1 ;
        k = k + numRegulatorsInFig ;
    end
end

output = struct;

output.regulatorNames = regulatorNames ;
output.nregulators = nregulators ;
output.targetName = targetName ;
output.experiments = datasets ;
output.nexps = nexps ;
output.nSplineSamples = nBootSamples ;
output.hSteps = hSteps ;

output.timepointsAll = timepointsAll ;
output.timeAll = timeAll ;
output.timepointsColAll = timepointsColAll ;
output.regDataAll = regDataAll ;
output.regDataColAll = regDataColAll;
% output.smoothRegProfileAllBoot = smoothRegProfileAllBoot ;
output.smoothRegProfileAll = smoothRegProfileAll ;

output.targetDataAll = targetDataAll ;
output.expNameAll = expNameAll ;
output.ntimepointsAll = ntimepointsAll ;
output.nreplicatesAll = nreplicatesAll ;
output.regIDAll = regIDAll ;
output.targetIDAll = targetIDAll ;
output.eHatMatAll = eHatMatAll;
output.sampleSizeAll = sampleSizeAll ;

output.regProfileAll = regProfileAll ;
output.regProfileGradAll = regProfileGradAll ;
output.regProfileAbsGradAll = regProfileAbsGradAll ;
output.regProfileAllBoot = regProfileAllBoot;
% output.regProfileGradAllBoot = regProfileGradAllBoot;
output.regProfileAbsGradAllBoot = regProfileAbsGradAllBoot;
output.regProfileActivationAll= regProfileActivationAll ;
output.regProfileActivationAllProposed = regProfileActivationAllProposed ;
output.minSwTime= minSwTime ;

output.YAll = YAll ;
output.maxRegProfiles = maxRegProfiles ;
output.minRegProfiles = minRegProfiles ;
output.maxRegProfileBoot = maxRegProfileBoot ;
output.minRegProfileBoot = minRegProfileBoot ;
output.minTargetData = minTargetData ;
output.maxTargetData = maxTargetData ;
output.minData = minData ;
output.maxData = maxData ;
output.minDataReg = minDataReg ;
output.maxDataReg = maxDataReg ;

output.timeFuncSize = timeFuncSize ;

output.smoothRegProfileColAll = smoothRegProfileColAll ;
output.resAll = resAll ;

output.targetSmoothData = targetSmoothData ;

end