function output = createThresholdsPrior(profilesOutput,timeFuncSize,profileFuncSize,lBoundStandGrads)

experiments = profilesOutput.experiments ;
regulatorNames = profilesOutput.regulatorNames ;

nexps = size(experiments,1);
nregs = length(regulatorNames);

% nregsAll = length(splinesOutput.regulatorNames) ;

regProfileAll = profilesOutput.regProfileAll ;
% regProfileAbsGradAll = profilesOutput.regProfileAbsGradAll ;
regProfileAbsGradAllBoot = profilesOutput.regProfileAbsGradAllBoot ;

minRegProfiles = profilesOutput.minRegProfiles ;
maxRegProfiles = profilesOutput.maxRegProfiles ;

hSteps = profilesOutput.hSteps ; 

muregProfileAbsGradAll = zeros(nexps,nregs,timeFuncSize) ;
stdregProfileAbsGradAll = zeros(nexps,nregs,timeFuncSize) ;
coeffregProfileAbsGradAll = zeros(nexps,nregs,timeFuncSize) ;
for x=1:nexps
    for z=1:nregs
        for i=1:timeFuncSize
            regProfileAbsGradAllBootiz = squeeze(regProfileAbsGradAllBoot{x}(i,:,z)) ;
            muregProfileAbsGradAll(x,i,z) = mean(regProfileAbsGradAllBootiz) ;
            stdregProfileAbsGradAll(x,i,z) = std(regProfileAbsGradAllBootiz) ;
            coeffregProfileAbsGradAll(x,i,z) = muregProfileAbsGradAll(x,i,z)/stdregProfileAbsGradAll(x,i,z) ;
        end
    end
end
mucoeff = nanmean(coeffregProfileAbsGradAll(:)) ;
stdcoeff = nanstd(coeffregProfileAbsGradAll(:)) ;

mincoeffs = zeros(profileFuncSize-1,nregs) ;
coeffsSigmoid = zeros(profileFuncSize-1,nregs) ;
threshPriorDensity = zeros(profileFuncSize-1,nregs) ;

regProfileActivationXZ = zeros(timeFuncSize,1) ;
for z = 1:nregs
    le = linspace(minRegProfiles(z),maxRegProfiles(z),profileFuncSize) ;
    for i=1:profileFuncSize-1
        lei = le(i+1) ;
        coeffsZI = [] ;
        for x = 1:nexps
            regProfXZ = regProfileAll{x}(:,z) ;
             
            %%% get activation
            regProfileActivationXZ(hSteps(x):end) = getRegProfileActivation(regProfXZ(hSteps(x):end),lei);
            regProfileActivationXZ(1:hSteps(x)-1) = regProfileActivationXZ(hSteps(x)) ;
            
            %%% get vector with state at each (dense) time point for all parents
            twosMat = repmat(pow2(0:(1-1)),timeFuncSize,1);
            stateChanges = getGeneActivation(regProfileActivationXZ,twosMat);
            
            swPos = getSwitchPointsThresh(stateChanges) ;% states used and at which times
            if ~isempty(swPos)
                coeffsZI = [coeffsZI coeffregProfileAbsGradAll(x,swPos,z)]    ;
            end
        end
        if ~isempty(coeffsZI)
            mincoeffs(i,z) = min(coeffsZI) ;
            coeffsSigmoid(i,z) = normcdf(mincoeffs(i,z),mucoeff,stdcoeff)  ;  %
        end
    end
    gaps = find(coeffsSigmoid(:,z) == 0) ; % get the gap regions
    numGapRegions = length(find(diff(gaps) >1))+1 ;  % get the diff within gaps. If any is longer than one we have more than one gaps. The length +1 gives the number of gaps
    if isempty(gaps)         
        threshPriorDensity(:,z) = coeffsSigmoid(:,z)/trapz(le(2:end),(coeffsSigmoid(:,z))) ;
    elseif (numGapRegions == 1)
        w = [(le(min(gaps)) - le(1))/(le(end) - le(1)) ; (le(max(gaps)) - le(min(gaps)))/(le(end) - le(1)) ; (le(end) - le(max(gaps)))/(le(end) - le(1))] ;
        k = 1 ;
        w2 = w.^k/(sum(w.^k)) ;
        threshPriorDensity(1:min(gaps)-1,z) = w2(1)*coeffsSigmoid(1:min(gaps)-1,z)/trapz(le(2:min(gaps)),coeffsSigmoid(1:min(gaps)-1,z)) ;
        threshPriorDensity(min(gaps):max(gaps),z) = w2(2)/(le(max(gaps)+1) - le(min(gaps)+1)) ;
        threshPriorDensity(max(gaps)+1:end,z) = w2(3)*coeffsSigmoid(max(gaps)+1:end,z)/trapz(le(max(gaps)+2:end),coeffsSigmoid(max(gaps)+1:end,z)) ;
    else
        gapsFinishIdx = [find(diff(gaps) >1 ) ; length(gaps)] ; 
        gapsFinish = gaps(gapsFinishIdx) ;
        gapsStartIdx = [1 ; gapsFinishIdx(1:end-1)+1] ;
        gapsStart = gaps(gapsStartIdx) ;
        threshPriorDensity(1:gapsStart(1)-1,z) = ((le(gapsStart(1)) - le(1))/(le(end) - le(1)))*(coeffsSigmoid(1:gapsStart(1)-1,z)/trapz(le(2:gapsStart(1)),coeffsSigmoid(1:gapsStart(1)-1,z))) ;
        for i=1:length(gapsStart)-1
            threshPriorDensity(gapsStart(i):gapsFinish(i),z) = ((le(gapsFinish(i)) - le(gapsStart(i)))/(le(end) - le(1)))/(le(gapsFinish(i)+1) - le(gapsStart(i)+1)) ;
            threshPriorDensity(gapsFinish(i)+1:gapsStart(i+1)-1,z) = ((le(gapsStart(i+1)) - le(gapsFinish(i)+1))/(le(end) - le(1)))*coeffsSigmoid(gapsFinish(i)+1:gapsStart(i+1)-1,z)/trapz(le(gapsFinish(i)+2:gapsStart(i+1)),coeffsSigmoid(gapsFinish(i)+1:gapsStart(i+1)-1,z)) ;
        end
        threshPriorDensity(gapsStart(end):gapsFinish(end),z) = ((le(gapsFinish(end)) - le(gapsStart(end)))/(le(end) - le(1)))*(1/length(le(gapsStart(end):gapsFinish(end)))) ;
        threshPriorDensity(gapsFinish(end)+1:end,z) = ((le(end) - le(gapsFinish(end)))/(le(end) - le(1)))*coeffsSigmoid(gapsFinish(end)+1:end,z)/trapz(le(gapsFinish(end)+2:end),coeffsSigmoid(gapsFinish(end)+1:end,z)) ;
    end
end
threshPriorDensity = [zeros(1,nregs) ; threshPriorDensity] ;
mincoeffs = [zeros(1,nregs) ; mincoeffs] ;

if (strcmp(lBoundStandGrads{1},'fixed'))
    lBStandGrads = lBoundStandGrads{2} ;
    lBStandGrads = repmat(lBStandGrads,nregs,1) ;
elseif (strcmp(lBoundStandGrads{1},'dataPrctile'))
    mm = mincoeffs(:) ;
    mm = mm(mm~=0) ;
    lBStandGrads = prctile(mm,lBoundStandGrads{2}) ;
    lBStandGrads = repmat(lBStandGrads,nregs,1) ;
elseif (strcmp(lBoundStandGrads{1},'dataPrctileRegs'))
    lBStandGrads = zeros(nregs,1) ;
    for z=1:nregs
        mm = mincoeffs(:,z) ; 
        mm = mm(mm~=0) ;
        lBStandGrads(z) = prctile(mm,lBoundStandGrads{2}) ;
    end
end

threshPriorDensity2 = zeros(profileFuncSize,nregs) ;
for z = 1:nregs
    le = linspace(minRegProfiles(z),maxRegProfiles(z),profileFuncSize) ;
    for i=1:profileFuncSize
        if (mincoeffs(i,z) < lBStandGrads(z) && mincoeffs(i,z)~=0)            
            threshPriorDensity2(i,z) = 0 ;            
        else
            threshPriorDensity2(i,z) = threshPriorDensity(i,z) ;
        end
    end
    if sum(threshPriorDensity2(:,z))~=0
        threshPriorDensity2(:,z) = threshPriorDensity2(:,z)/trapz(le,threshPriorDensity2(:,z)) ;
    end
end

output.threshPriorDensity = threshPriorDensity ;
output.threshPriorDensity2 = threshPriorDensity2 ;

end


function [states] = getGeneActivation(regExpr,twos)

states = binaryMatToRow(regExpr,twos);

end

function [rows] = binaryMatToRow(binMatrix,twos)

rows = sum(twos .* binMatrix,2)+1;

end

function swPos = getSwitchPointsThresh(activation)

swPos = find(abs(activation(2:end) - activation(1:(end-1))))+1;

end

function [activation] = getRegProfileActivation(expression,threshold)

if (isempty(expression) == 0)
    activation = expression >= threshold;
else
    [l,~] = size(expression) ;
    activation = zeros(l,1) ;
end

end


