function output = createProfiles_Priors_Proposals(data_filenames,basenet_filename,smoothProfiles_filename,targetName,profiles,priors,nregsPriors,regsPriors,threshPriors,varargin)

currentpath = cd('..');
parentpath = pwd();
cd(currentpath);

datasets = cell(size(data_filenames,2),1);
for i=1:size(data_filenames,2)
    load([parentpath '/output/' data_filenames{i}],'dataset') ;  
    datasets{i} = dataset ;
end

load([parentpath '/output/' basenet_filename],'net') ;

regNames = net.netForInterType{strcmp(net.netForInterType{:,3},targetName),1} ;

load([parentpath '/output/' smoothProfiles_filename],'smoothProfiles')

%% default variables
maxNumRegs = [] ; % max number of regulators allowed by the prior distributions
hyperparamsLambda0 = 0.125 ; % % 2; % 0.125 ; mean for poisson numRegulators Prior

regulatorsPriorType = 'range' ; % type of regulators prior

thresholdsPriorType = 'gradient' ; % type of thresholds prior
thresholdsPriorLBoundStandGrads = {'dataPrctile',25} ; % lower cut off point for the gradient computed as a percentile of the gradient values across candidate regulators

hyperparamsDeg = 0.345 ;  % mean of the prior for regulator degradation
hyperparamsMinDegRate = 0.1; % minimum value for the truncated proposal ratio for deg
hyperparamsMaxDegRate =  2; % maximum value for the truncated proposal ratio for deg
hyperparamsPrior_deg_nu0 = 10 ; % degrees of freedom for the non-central chi2 prior distribution for deg
hyperparamsPrior_deg_s02 = 1/0.345 ; % scale parameter for the non-central chi2 prior distribution for deg
hyperparamsPrecDeg0 = 10^(-3) ; % degrees of freedom for the chi2 prior distribution for precision
hyperparamsPrecSigma0 = 10^(-3) ;  % scale for the chi squared prior distribution for precision
profileVecSize = 1000 ;  % vector size for the dense partition of regulators profile region
priorPlotsOn = 0 ; % (0/1) for plotting priors' figure
regProfilePlotsOn = 0 ; % (0/1) for plotting regulator profiles figure
numRegulatorsInFig = 4 ; % num of regulators to be plotted in each figure
minSwTimeDefault = 1 ; % (0/1) default sets this equal to delaytime
minSwTime = [] ; % if 0, set a fixed vector of length numExps 

%% processing optional inputs
if mod(length(varargin),2)
    error('Syntax error on command line. Enter Name Value pairs');
end
for i=1:2:length(varargin)-1
    switch varargin{i}
        case {'maxNumRegs'}
            maxNumRegs = varargin{i+1};
        case('regulatorsPriorType')
            regulatorsPriorType = varargin{i+1};
        case('thresholdsPriorType')
            thresholdsPriorType = varargin{i+1};
        case('thresholdsPriorLBoundStandGrads')
            thresholdsPriorLBoundStandGrads = varargin{i+1};
        case('hyperparamsDeg')
            hyperparamsDeg = varargin{i+1};
        case('hyperparamsMinDegRate')
            hyperparamsMinDegRate = varargin{i+1};
        case('hyperparamsMaxDegRate')
            hyperparamsMaxDegRate = varargin{i+1};
        case('hyperparamsPrior_deg_nu0')
            hyperparamsPrior_deg_nu0 = varargin{i+1};
        case('hyperparamsPrior_deg_s02')
            hyperparamsPrior_deg_s02 = varargin{i+1};
        case('hyperparamsLambda0')
            hyperparamsLambda0 = varargin{i+1};
        case('hyperparamsPrecDeg0')
            hyperparamsPrecDeg0 = varargin{i+1};
        case('hyperparamsPrecSigma0')
            hyperparamsPrecSigma0 = varargin{i+1};  
        case('profileVecSize')
            profileVecSize = varargin{i+1};  
        case('priorPlotsOn')
            priorPlotsOn = varargin{i+1}; 
        case('regProfilePlotsOn')
            regProfilePlotsOn = varargin{i+1}; 
        case('numRegulatorsInFig')
            numRegulatorsInFig = varargin{i+1};
        case('minSwTimeDefault')
            minSwTimeDefault = varargin{i+1}; 
        case('minSwTime')
            minSwTime = varargin{i+1}; 
        otherwise
            error(['Unknown property ' varargin{i}  ' in the command line']);
    end
end

% if isempty(profiles)
    profiles = getProfiles(targetName,regNames,datasets,smoothProfiles,regProfilePlotsOn,numRegulatorsInFig,minSwTimeDefault,minSwTime); 
% end

tic; priorsPropOutput = createPriorsProp(profiles,priors,nregsPriors,regsPriors,threshPriors,...
    maxNumRegs,hyperparamsLambda0,...
    regulatorsPriorType,thresholdsPriorType,thresholdsPriorLBoundStandGrads,...
    hyperparamsDeg,hyperparamsMinDegRate,hyperparamsMaxDegRate,hyperparamsPrior_deg_nu0,...
    hyperparamsPrior_deg_s02,hyperparamsPrecDeg0,hyperparamsPrecSigma0,...
    profileVecSize,priorPlotsOn,numRegulatorsInFig) ; toc ;

output = struct ;
output.profilesOutput = profiles ;
output.priorsPropOutput = priorsPropOutput ;

end