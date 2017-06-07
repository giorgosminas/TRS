function output = createSmoothProfiles(basenet_filename,data_filenames,varargin)

% create the smooth profiles of the regulators
% optionally creates Wild Bootstrap samples of the profiles (Wu, 1986)
% It has a set of optional input variables that can be set by the 
% name,value format (see section default parameters and processing optional input)

currentpath = cd('..');
parentpath = pwd();
cd(currentpath);

rng('shuffle')

dataset = cell(size(data_filenames,2),1) ;
for i=1:size(data_filenames,2)
    dataset{i} = load([parentpath '/output/' data_filenames{i}],'dataset') ;    
end

load([parentpath '/output/' basenet_filename],'net') ;

%% default parameters
regulatorType = 0 ; % 0: delayedMRNA, 1: protein interp(olated from mRNA), 2: protein observ(ations)   
nBootSamples = 100 ; % # of spline samples to be generated for bootstrap
plotsOn = 1 ;  % 1: plots to be generated
numParentsInFig = 6 ;  % # of parent profiles in each fig
smoothPar = [0.1 0.1] ; % smoothing  parameters, see MATLAB function fit, option 'smoothingspline'
delayTime = 1 ;  % if profiles are expression data and regulator=0, the time of delay applied to expression data
deg_p = log(2)/2 ; % if profiles are expression data and regulator=1, the degradation parameter for proteins
trans_p = log(2)/1 ; % if profiles are expression data and regulator=1, the translation parameter
funcSize = 100; % # of timepoints of the smoothed parameter

%% processing optional input
if mod(length(varargin),2)
    error('Syntax error on command line. Enter Name Value pairs');
end
for i=1:2:length(varargin)-1
    switch varargin{i}
        case ('regulator')
            regulatorType = varargin{i+1};
        case ('nBootSamples')
            nBootSamples = varargin{i+1};
        case ('plotsOn')
            plotsOn = varargin{i+1};
        case ('numParentsInFig')
            numParentsInFig = varargin{i+1};
        case ('smoothPar')
            smoothPar = varargin{i+1};
        case ('delayTime')
            delayTime = varargin{i+1};
        case ('deg_p')            
            deg_p = varargin{i+1};
        case ('trans_p')            
            trans_p = varargin{i+1};
        case ('funcSize')            
            funcSize = varargin{i+1};
        otherwise
            error(['Unknown property ' varargin{i}  ' in the command line']);
    end
end


if ischar(regulatorType)
    if strcmp(regulatorType,'delayedMRNA')
        regulatorType = 0 ;
    elseif strcmp(regulatorType,'protein interp')
        regulatorType = 1 ;
    elseif strcmp(regulatorType,'protein observations')
        regulatorType = 2 ;
    end
elseif ~ismember(regulatorType,[0 1 2])
    error('regulator not specified correctly');
end

if ~isempty(smoothPar)
    if ~isnumeric(smoothPar)
        error('smoothPar should be numeric');
    elseif ~ismember(sum(size(smoothPar)),[2 3])
        error('smoothPar should be either scalar or a 2x1 vector');
    elseif isscalar(smoothPar)
        if smoothPar>1 || smoothPar<0
            error('smoothPar should be between 0 and 1');
        else
            smoothPar = smoothPar*ones(1,2) ;
        end
    elseif any(smoothPar>1) ||  any(smoothPar<0)
        error('smoothPar should be between 0 and 1');
    end
end

if funcSize<=0
    error('funcSize should be a positive integer');    
else
   funcSize = ceil(funcSize);
end

if ~isnumeric(nBootSamples) || nBootSamples<1
    error('nSplineSamples should be numeric >=1');
elseif ~ismember(plotsOn,[0 1])
    error('plotsOn should be 0 or 1');
elseif ~isnumeric(numParentsInFig) || numParentsInFig<1
    error('numParentsInFig should be numeric >=1');
elseif ~isnumeric(delayTime) || delayTime<0
    error('delayTime should be numeric >=0');
elseif ~isnumeric(deg_p) || deg_p<=0
    error('deg_p should be numeric >0');
elseif ~isnumeric(trans_p) || trans_p<=0
    error('trans_p should be numeric >0');
end

%% get various variables
regulators = net.regulatorNamesForInterType ;
regulatorNames = regulators(:,1) ;

nexps = length(dataset);
nregs = length(regulatorNames);

%% check to see if we have all of the regulators in all experiments
hasData = zeros(nregs,nexps);
for x = 1:nregs
    p = regulatorNames{x};
    for y = 1:nexps
        e = dataset{y}.dataset;
        hasData(x,y) = findIndex(p,e.orderedNames);
    end
end
foundData = sum(hasData > 0,2) == nexps;
if sum(foundData) < nregs
    disp('Removing missing regulators');
    regulatorNames = regulatorNames(foundData);
    nregs = length(regulatorNames);
    par_str = regulatorNames{1};
    for x = 2:length(regulatorNames)
        par_str = [par_str ', ' regulatorNames{x}]; %#ok<AGROW>
    end
    if length(regulatorNames) == 1
        disp(['Parent probe (1): ' par_str]);
    else
        disp(['Parent probes (' int2str(length(regulatorNames)) '): ' par_str]);
    end
end


%% create space for the experiments' variables
expNameAll = cell(nexps,1); % name of each experiment
ntimepointsAll = zeros(nexps,1); % # of timepoints in each experiment
replicatesAll = zeros(nexps,1); % # of replicates in each experiment
timepointsAll = cell(nexps,1); % timepoints in each experiment
timeAll = cell(nexps,1);  % times for interpolation of the smooth profile
timepointsColAll = cell(nexps,1) ;  % column of replicated timepoints
regDataAll = cell(nexps,nregs);  % regulator data
regDataColAll = cell(nexps,nregs);  % regulator data in column
smoothRegProfileAllBoot = cell(nexps,nregs); % Bootstrap samples of smoothed observed (mRNA or protein) profile of the regulator
smoothRegProfileAll = cell(nexps,nregs); % (median) of the above
smoothRegProfileColAll = cell(nexps,nregs) ; % as above but in column
regIDAll = zeros(nexps,nregs); % idx of row containing regulator profile in data

%% and for regulators expression, activation, etc.
regProfileAll = cell(nexps,1) ; % profiles of regulators
regProfileGradAll = cell(nexps,1) ; % their gradient
regProfileAbsGradAll = cell(nexps,1) ; % the absolute value of their gradient 
regProfileAllBoot = cell(nexps,1);  % bootstrap samples of reg profiles
regProfileGradAllBoot = cell(nexps,1);  % their gradient 
regProfileAbsGradAllBoot = cell(nexps,1); % their absolute gradient 

regProfileActivationAll = cell(nexps,1);  % activation functions of regulators
regProfileActivationAllProposed = cell(nexps,1);  % activation functions of proposed regulators (to be used in mcmc)

%% and min and max of various variables
minRegProfilesBoot = 10^6*ones(nBootSamples,nregs) ; % min smooth reg profile over all bootstrap samples
maxRegProfilesBoot = -10^6*ones(nBootSamples,nregs) ; % max

minRegProfiles = 10^20*ones(1,nregs) ; % min of reg profiles
maxRegProfiles = -10^20*ones(1,nregs) ; % max

minDataReg = 10^20*ones(1,nregs) ; % min of data
maxDataReg = -10^20*ones(1,nregs) ; % max 

minSmoothRegProfileBoot = [] ; % min smooth reg profile over all bootstrap samples
maxSmoothRegProfileBoot = [] ; % max

resAll = cell(nexps,nregs) ; % residuals in wild bootstrap
hSteps = zeros(nexps,1) ; % idx of first point in regulator profileafter delay (for delayed mRNA)

%% get variables
for x = 1:nexps
    %% get values of experiment variables
    experiment = dataset{x}.dataset;
    
    expName = experiment.name;
    ntimepoints = experiment.ntimepoints; 
    nreplicates = experiment.nreplicates; 
    
    expNameAll{x} = expName;
    ntimepointsAll(x) = ntimepoints;
    replicatesAll(x) = nreplicates;
        
    timepoints = experiment.timepoints;
    timepoints = timepoints - timepoints(1); 
    timepointsAll{x} = timepoints;
    
    timepointsCol = repmat(timepoints,1,nreplicates)' ;
    timepointsColAll{x} = timepointsCol ;    
    
    time = linspace(0,timepoints(end),funcSize);
    timeAll{x} = time;
    
    regProfileAll{x} = nan(funcSize,nregs) ;
    regProfileGradAll{x} = nan(funcSize,nregs) ;
    regProfileAbsGradAll{x} = nan(funcSize,nregs) ;
    
    regProfileAllBoot{x} = nan(funcSize,nBootSamples, nregs);
    regProfileGradAllBoot{x} = nan(funcSize,nBootSamples, nregs);
    regProfileAbsGradAllBoot{x} = nan(funcSize,nBootSamples, nregs);
    
    regProfileActivationAll{x} = nan(funcSize,nregs);
    regProfileActivationAllProposed{x} = nan(funcSize,nregs);

    regProfileX = regProfileAll{x} ;
    regProfileGradX = regProfileGradAll{x} ;
    regProfileAbsGradX = regProfileAbsGradAll{x} ;
    
    regProfileXBoot = regProfileAllBoot{x};
    regProfileGradXBoot = regProfileGradAllBoot{x} ;
    regProfileAbsGradXBoot = regProfileAbsGradAllBoot{x} ;
    for z = 1:nregs
        
        %%% get the parent ID to get its data
        regID = findIndex(regulatorNames{z},experiment.orderedNames);
        regIDAll(x,z) = regID;
        
        dataCol = experiment.data(regID,:)' ;
        regDataColAll{x,z} = dataCol;
        
        minDataReg(z) = min([minDataReg(z) min(dataCol)]) ;
        maxDataReg(z) = max([maxDataReg(z) max(dataCol)]) ;
        
        d = reshape(dataCol',ntimepoints,nreplicates);
        data = d';
        regDataAll{x,z} = data;
        
        %%%% fit a smoothed spline for continuous mRNA
        if ~isempty(smoothPar)
            [f,~,out] = fit(timepointsCol,dataCol,'smoothingspline','SmoothingParam', smoothPar(1));
        else
            [f,~,out] = fit(timepointsCol,dataCol,'smoothingspline');
        end
        smoothRegProfile = zeros(funcSize,nBootSamples) ;
        smoothRegProfile(:,1) = feval(f,time);
        smoothRegProfileAll{x,z} = smoothRegProfile ;
            
        if nBootSamples>1
            
            smoothRegProfileCol = feval(f,timepointsCol);
            smoothRegProfileColAll{x,z} = smoothRegProfileCol;
            res = out.residuals ;
            resAll{x,z} = res ;
            
            for n=2:nBootSamples
                sampledData = zeros(length(res),1) ;
                for q = 1:nreplicates
                    sampledRes = res((q-1)*ntimepoints+1:q*ntimepoints).*randn(ntimepoints,1) ; % randsample(res((q-1)*timepoints+1:q*timepoints),timepoints,true) ;
                    sampledData((q-1)*ntimepoints+1:q*ntimepoints) = smoothRegProfileCol((q-1)*ntimepoints+1:q*ntimepoints) + sampledRes;
                end
                sampledData = max(0,sampledData) ;  
                if ~isempty(smoothPar)
                    f = fit(timepointsCol,sampledData,'smoothingspline', 'SmoothingParam', smoothPar(2));
                else
                    f = fit(timepointsCol,dataCol,'smoothingspline');
                end
                smoothRegProfile(:,n) = feval(f,time);
            end
            
            smoothRegProfileAllBoot{x,z} = smoothRegProfile ;
            smoothRegProfileAll{x,z} = median(smoothRegProfile,2) ;
            
            minSmoothRegProfileBoot = min([min(min(smoothRegProfile)) minSmoothRegProfileBoot]) ;
            maxSmoothRegProfileBoot = max([max(max(smoothRegProfile)) maxSmoothRegProfileBoot]) ;
        end
        
        %%% if regulator is protein get regProfile via protein ODE on
        %%% smoothed profile
        if regulatorType==1
                        
            regProfileX(1,z) = smoothRegProfileAll{x,z}(1) ;
            for q=2:funcSize
                xq = exp(-deg_p*time(q)) ;
                x1q = exp(deg_p*time(1:q)) ;
                regProfileX(q,z) = regProfileX(1,z)*xq + trans_p*xq*myTrapz(time(1:q),smoothRegProfileAll{x,z}(1:q).*x1q') ;
            end
            
            minRegProfiles = min([squeeze(min(regProfileX)) ; minRegProfiles]) ;
            maxRegProfiles = max([squeeze(max(regProfileX)) ; maxRegProfiles]) ;
            
            p = regProfileX(:,z) ;
            regProfileGradX(:,z) = gradient(p,time(2))';
            regProfileAbsGradX(:,z) = abs(regProfileGradX(:,z));
            
            hSteps(x) = 1 ;
            for n=1:nBootSamples
                
                p = nan(funcSize,1);
                p(1) = smoothRegProfile(1,n) ;
                for q=2:funcSize
                    xq = exp(-deg_p*time(q)) ;
                    x1q = exp(deg_p*time(1:q)) ;
                    p(q) = p(1)*xq + trans_p*xq*myTrapz(time(1:q),smoothRegProfile(1:q,n).*x1q') ;
                end
                
                regProfileXBoot(hSteps(x):end,n,z) = p ;
                
                minRegProfilesBoot(n,:) = min([squeeze(min(regProfileXBoot(:,n,:)))' ; minRegProfilesBoot(n,:)]) ;
                maxRegProfilesBoot(n,:) = max([squeeze(max(regProfileXBoot(:,n,:)))' ; maxRegProfilesBoot(n,:)]) ;
                
                regProfileGradXBoot(:,n,z) = gradient(p,time(2))';
                regProfileAbsGradXBoot(:,n,z) = abs(regProfileGradXBoot(:,n,z));
            end
            
            
        elseif ismember(regulatorType,[0 2])
            
            if regulatorType==2 
                delayTime = 0 ;
                hSteps(x) = 1 ;
            else
                hSteps(x) = find(time<=delayTime,1,'last') ;
            end
            
            regProfileX(hSteps(x):end,z) = smoothRegProfileAll{x,z}(1:funcSize - hSteps(x)+1) ;
            
            minRegProfiles = min([squeeze(min(regProfileX)) ; minRegProfiles]) ;
            maxRegProfiles = max([squeeze(max(regProfileX)) ; maxRegProfiles]) ;
            
            p = regProfileX(:,z) ;
            regProfileGradX(:,z) = gradient(p,time(2))';
            regProfileAbsGradX(:,z) = abs(regProfileGradX(:,z));
            
            for n=1:nBootSamples
                regProfileXBoot(hSteps(x):end,n,z) = smoothRegProfile(1:funcSize - hSteps(x)+1,n) ;
                
                minRegProfilesBoot(n,:) = min([squeeze(min(regProfileXBoot(:,n,:)))' ; minRegProfilesBoot(n,:)]) ;
                maxRegProfilesBoot(n,:) = max([squeeze(max(regProfileXBoot(:,n,:)))' ; maxRegProfilesBoot(n,:)]) ;
                
                p = regProfileXBoot(:,n,z) ;
                regProfileGradXBoot(:,n,z) = gradient(p,time(2))';
                regProfileAbsGradXBoot(:,n,z) = abs(regProfileGradXBoot(:,n,z));
            end
        end
    end
    regProfileAll{x} = regProfileX ;
    regProfileGradAll{x} = regProfileGradX ;
    regProfileAbsGradAll{x} = regProfileAbsGradX ;
    
    regProfileAllBoot{x} = regProfileXBoot ;
    regProfileGradAllBoot{x} = regProfileGradXBoot ;
    regProfileAbsGradAllBoot{x} = regProfileAbsGradXBoot ;
end
minData = min(minDataReg) ;
maxData = max(maxDataReg) ;

if plotsOn
    counter = nregs ;
    k = 0 ;
    while (counter > 0)
        m1 = min(nregs-k,numParentsInFig) ;
        figure('name','Regulators data','numbertitle','off','windowstyle','docked');
        lmargin = 0.125 ;
        bmarginw = 0.04 ;
        rmargin = 0.05 ;
        plotWidth = (1 - lmargin - nexps*bmarginw - rmargin)/nexps ;
        
        umargin = 0.05 ;
        bmarginh = 0.05 ;
        dmargin = 0.05 ;
        plotHeight = (1 - umargin - (m1-1)*bmarginh - dmargin)/m1 ;
        
        for x=1:nexps
            for z=k+1:k+m1
                positionVector = [lmargin + (x-1)*(bmarginw + plotWidth) dmargin + (m1-(z-k))*(bmarginh + plotHeight)  plotWidth  plotHeight] ;
                subplot('Position',positionVector)
                plot(repmat(timepointsAll{x},1,replicatesAll(x)),regDataColAll{x,z},'ob');
                hold all
                plot(timeAll{x},smoothRegProfileAllBoot{x,z}(:,1),'-k')                
                if (z == 1)
                    title(['Experiment ' num2str(x)],'FontSize',12.5,'FontWeight','bold')
                end
                if (x == 1)
                    ylabel(['reg ' num2str(z)],'FontSize',12,'FontWeight','bold')
                end
                
                for n=2:nBootSamples
                    plot(timeAll{x},smoothRegProfileAllBoot{x,z}(:,n),'-c')
                end
                plot(timeAll{x},smoothRegProfileAll{x,z},'-.k','lineWidth',1.5)
                axis tight
            end
        end
        counter = counter - m1 ;
        k = k + numParentsInFig ;
    end
end

output = struct;

output.regulatorNames = regulatorNames ;
output.nregulators = nregs ;
output.nexps = nexps ;
output.regulator = regulatorType ;
output.nBootSamples = nBootSamples ;
output.regDataAll = regDataAll ;
output.regDataColAll = regDataColAll;
output.regIDAll = regIDAll ;
output.regProfileAll = regProfileAll ;
output.regProfileGradAll = regProfileGradAll ;
output.regProfileAbsGradAll = regProfileAbsGradAll ;
output.regProfileAllBoot = regProfileAllBoot;
output.regProfileGradAllBoot = regProfileGradAllBoot;
output.regProfileAbsGradAllBoot = regProfileAbsGradAllBoot;
output.regProfileActivationAll= regProfileActivationAll ;
output.regProfileActivationAllProposed = regProfileActivationAllProposed ;
output.hSteps = hSteps ;
output.maxRegProfiles = maxRegProfiles ;
output.minRegProfiles = minRegProfiles ;
output.minRegProfilesBoot = minRegProfilesBoot ;
output.maxRegProfilesBoot = maxRegProfilesBoot ;
output.minData = minData ;
output.maxData = maxData ;
output.minDataReg = minDataReg ;
output.maxDataReg = maxDataReg ;
output.funcSize = funcSize ;
output.smoothRegProfileAllBoot = smoothRegProfileAllBoot ;
output.smoothRegProfileAll = smoothRegProfileAll ;
output.smoothRegProfileColAll = smoothRegProfileColAll ;
output.resAll = resAll ;

end



function z = myTrapz(x,y,dim)
%TRAPZ  Trapezoidal numerical integration.
%   Z = TRAPZ(Y) computes an approximation of the integral of Y via
%   the trapezoidal method (with unit spacing).  To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, TRAPZ(Y) is the integral of Y. For matrices, TRAPZ(Y)
%   is a row vector with the integral over each column. For N-D
%   arrays, TRAPZ(Y) works across the first non-singleton dimension.
%
%   Z = TRAPZ(X,Y) computes the integral of Y with respect to X using
%   the trapezoidal method.  X and Y must be vectors of the same
%   length, or X must be a column vector and Y an array whose first
%   non-singleton dimension is length(X).  TRAPZ operates along this
%   dimension.
%
%   Z = TRAPZ(X,Y,DIM) or TRAPZ(Y,DIM) integrates across dimension DIM
%   of Y. The length of X must be the same as size(Y,DIM)).
%
%   Example: If Y = [0 1 2
%                    3 4 5]
%
%   then trapz(Y,1) is [1.5 2.5 3.5] and trapz(Y,2) is [2
%                                                       8];
%
%   Class support for inputs X, Y:
%      float: double, single
%
%   See also SUM, CUMSUM, CUMTRAPZ, QUAD.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.16.4.6 $  $Date: 2006/12/15 19:27:35 $

%   Make sure x and y are column vectors, or y is a matrix.

perm = []; nshifts = 0;
if nargin == 3 % trapz(x,y,dim)
    perm = [dim:max(ndims(y),dim) 1:dim-1];
    y = permute(y,perm);
    m = size(y,1);
elseif nargin==2 && isscalar(y) % trapz(y,dim)
    dim = y; y = x;
    perm = [dim:max(ndims(y),dim) 1:dim-1];
    y = permute(y,perm);
    m = size(y,1);
    x = 1:m;
else % trapz(y) or trapz(x,y)
    if nargin < 2, y = x; end
    %   size(y)
    %   [y,nshifts] = myShiftdim(y);
    %   size(y)
    %   nshifts
    nshifts = 0;
    m = size(y,1);
    if nargin < 2, x = 1:m; end
end
if ~isvector(x)
    error('MATLAB:trapz:xNotVector', 'X must be a vector.');
end
x = x(:);
if length(x) ~= m
    if isempty(perm) % dim argument not given
        error('MATLAB:trapz:LengthXmismatchY',...
            'LENGTH(X) must equal the length of the first non-singleton dimension of Y.');
    else
        error('MATLAB:trapz:LengthXmismatchY',...
            'LENGTH(X) must equal the length of the DIM''th dimension of Y.');
    end
end

% The output size for [] is a special case when DIM is not given.
if isempty(perm) && isequal(y,[])
    z = zeros(1,class(y));
    return;
end

%   Trapezoid sum computed with vector-matrix multiply.
z = diff(x,1,1).' * (y(1:m-1,:) + y(2:m,:))/2;

siz = size(y); siz(1) = 1;
z = reshape(z,[ones(1,nshifts),siz]);
if ~isempty(perm), z = ipermute(z,perm);
end
end





