% TRS, Transcriptional Regulation Switch 
% by 
% Copyright (C) 2016 Giorgos Minas, Dafyd Jenkins, David Rand, Barbel Finkenstadt
% and University of Warwick
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% clear all  %#ok<CLSCR>

% run from directory '.../trRegCode/scripts'
cd('..'); 
parentpath = pwd();
cd('functions/')

global dataset_filenames net_filename smoothProfiles_filename targetName priors_filename namestr
if isempty(dataset_filenames) || isempty(net_filename) || isempty(smoothProfiles_filename)
     dlg_title = 'Filenames and Target' ;
     prompt = {'data filename','data filename','data filename','data filename','data filename','data filename',...
         'Base network filename','smoothProfiles_filename','targetName'} ;
     num_lines = 1 ;
     defaultPreC_dlg = {'','','','','','','','',''} ;
     precomputePriorsdlg = inputdlg(prompt,dlg_title,num_lines,defaultPreC_dlg) ;
     dataset_filenames=[];
     for i=1:length(precomputePriorsdlg)-3
         if ~isempty(precomputePriorsdlg{i})
             dataset_filenames{i}=precomputePriorsdlg{i};
         end
     end
     net_filename = precomputePriorsdlg{end-2};
     smoothProfiles_filename = precomputePriorsdlg{end-1};
     targetName = precomputePriorsdlg{end};
elseif isempty(targetName)
     dlg_title = 'Target name';
     prompt = {'targetName'} ;
     num_lines = 1 ;
     defaultPreC_dlg = {''} ;
     precomputePriorsdlg = inputdlg(prompt,dlg_title,num_lines,defaultPreC_dlg) ;
     targetName = precomputePriorsdlg{1};
end

%% priors
profiles = [] ;
priors = [] ;
nregsPriors = 0 ;
regsPriors = 0 ;
threshPriors = 0 ;
degPrior = 0 ;
precPrior = 0 ;
prompt = {'use pre-computed prior','filename','number of regs prior as in file(0/1)',...
    'set of regs prior as in file(0/1)',...
    'threshold of regs prior as in file(0/1)'} ;
dlg_title = 'pre-computed priors' ;
num_lines = 1 ;
defaultPreC_dlg = {'0','','0','0','0'} ;
precomputePriorsdlg = inputdlg(prompt,dlg_title,num_lines,defaultPreC_dlg) ;
if str2double(precomputePriorsdlg{1})
    nregsPriors = str2double(precomputePriorsdlg{3}) ;
    regsPriors = str2double(precomputePriorsdlg{4}) ;
    threshPriors = str2double(precomputePriorsdlg{5}) ;
    priorsFileName = precomputePriorsdlg{2} ;
    load([parentpath '/output/' priorsFileName])
    profiles = profilesPriors.profilesOutput ;
    priors = profilesPriors.priorsPropOutput.priors ;
end
clear prompt dlg_title num_lines 
prompt = {'lambda parameter value for the poisson prior',...
          'max number of regulators',...
          'regulatorsPriorType (uniform or range)',...
          'thresholdsPriorType (uniform or gradient)',...
          'type of cut-off value for gradient (fixed, dataPrctile or dataPrctileRegs)',...
          'cut-off parameter',...
          'mean for the gamma prior of degradation rate',...
          'df for the chi2 prior distribution of precision',...
          'scale for the chi2 prior distribution of precision'
          };
dlg_title = 'priors' ;
num_lines=1;
if nregsPriors
    default_dlg{1} = num2str(hyperparamsLambda0) ;
    default_dlg{2} = num2str(maxNumRegs) ;
else
    default_dlg{1} = num2str(0.15) ;
    default_dlg{2} = num2str(8) ;
end
if regsPriors
    default_dlg{3} = regulatorsPriorType ;
else
    default_dlg{3} = 'range' ;
end
if threshPriors
    default_dlg{4} = thresholdsPriorType ;
    default_dlg{5} = thresholdsPriorLBoundStandGrads{1} ;
    default_dlg{6} = num2str(thresholdsPriorLBoundStandGrads{2}) ;
else
    default_dlg{4} = 'gradient' ;
    default_dlg{5} = 'dataPrctile' ;
    default_dlg{6} = '50' ;
end
default_dlg{7} = num2str(0.345) ; default_dlg{8} = num2str(1e-3) ; default_dlg{9} = num2str(1e-3) ;
priorsdlg = inputdlg(prompt,dlg_title,num_lines,default_dlg) ;
warningNumRegs = 0 ; warningRegs = 0 ; warningThresh = 0 ; 
if str2double(precomputePriorsdlg{1})
    if nregsPriors
        if str2num(priorsdlg{1})~=hyperparamsLambda0 %#ok<*ST2NM>
            hyperparamsLambda0 = str2num(priorsdlg{1}) ;
            warningNumRegs = 1 ;
        end
        if str2num(priorsdlg{2})~=maxNumRegs
            maxNumRegs = str2num(priorsdlg{2}) ;
            warningNumRegs = 1 ;
        end
    else
        hyperparamsLambda0 = str2num(priorsdlg{1}) ;
        maxNumRegs = str2num(priorsdlg{2}) ;
    end
    if regsPriors
        if ~strcmp(priorsdlg{3},regulatorsPriorType)
            regulatorsPriorType = priorsdlg{3} ;
            warningRegs = 1 ;
        end
    else
        regulatorsPriorType = priorsdlg{3} ;
    end
    if threshPriors
        if ~strcmp(priorsdlg{4},thresholdsPriorType)
            thresholdsPriorType = priorsdlg{4} ;
            warningThresh = 1 ;
        end
        if ~strcmp(priorsdlg{5},thresholdsPriorLBoundStandGrads{1})
            warningThresh = 1 ;
        end
        if str2num(priorsdlg{6})~=thresholdsPriorLBoundStandGrads{2} 
            thresholdsPriorLBoundStandGrads{2} = str2num(priorsdlg{6}) ;
            warningThresh = 1 ;
        end
    else
        thresholdsPriorType = priorsdlg{4} ;
        thresholdsPriorLBoundStandGrads{1} = priorsdlg{5} ;
        thresholdsPriorLBoundStandGrads{2} = str2num(priorsdlg{6}) ;
    end
else
    hyperparamsLambda0 = str2num(priorsdlg{1}) ;
    maxNumRegs = str2num(priorsdlg{2}) ;
    regulatorsPriorType = priorsdlg{3} ;
    thresholdsPriorType = priorsdlg{4} ;
    thresholdsPriorLBoundStandGrads{1} = priorsdlg{5} ;
    thresholdsPriorLBoundStandGrads{2} = str2num(priorsdlg{6}) ;
end

warning_str = 'The parameters of the priors for: ' ;
if warningNumRegs
    nregsPriors = 0 ;
    warning_str = [warning_str 'number of regs, '] ;    
end
if warningRegs
    regsPriors = 0 ;
    warning_str = [warning_str 'set of regs, '] ;    
end
if warningThresh
    threshPriors = 0 ;
    warning_str = [warning_str 'threshold of regs '] ;    
end
warning_str = [warning_str 'are changed and the corresponding priors will be re-computed'] ;
    
if any([warningNumRegs warningRegs warningThresh])
    warndlg(warning_str)
end

hyperparamsDeg = str2num(priorsdlg{7}) ;
hyperparamsPrecDeg0 = str2num(priorsdlg{8}) ;
hyperparamsPrecSigma0 = str2num(priorsdlg{9}) ;

hyperparamsMinDegRate = 0.1; % minimum value for the truncated proposal ratio for deg
hyperparamsMaxDegRate =  2; % maximum value for the truncated proposal ratio for deg
hyperparamsPrior_deg_nu0 = 10 ; % degrees of freedom for the non-central chi2 prior distribution for deg
hyperparamsPrior_deg_s02 = 1/0.345 ; % scale parameter for the non-central chi2 prior distribution for deg
clear warning_str precomputePriorsdlg prompt dlg_title num_lines warndlg warningNumRegs warningRegs warningThresh priorsdlg

%% other computational parameters
profileVecSize = 1000 ; % vector size for the dense partition of regulators profile region
priorPlotsOn = 1 ; % (0/1) for plotting priors' figure
regProfilePlotsOn = 0 ; % (0/1) for plotting regulator profiles figure
numRegulatorsInFig = 4 ; % num of regulators to be plotted in each figure
%%% min time for switches. earlier switches are cancelled out
minSwTimeDefault = 1 ; % (0/1) default sets this equal to delaytime
minSwTime = [] ; % if 0, set a fixed vector of length numExps 

profilesPriors = createProfiles_Priors_Proposals(dataset_filenames,net_filename,...
    smoothProfiles_filename,targetName,...
    profiles,priors,nregsPriors,regsPriors,threshPriors,...
    'maxNumRegs',maxNumRegs,...
    'hyperparamsLambda0',hyperparamsLambda0,...
    'regulatorsPriorType',regulatorsPriorType,...
    'thresholdsPriorType',thresholdsPriorType,...
    'thresholdsPriorLBoundStandGrads',thresholdsPriorLBoundStandGrads,...
    'hyperparamsDeg',hyperparamsDeg,'hyperparamsMinDegRate',hyperparamsMinDegRate,...
    'hyperparamsMaxDegRate',hyperparamsMaxDegRate,'hyperparamsPrior_deg_nu0',hyperparamsPrior_deg_nu0,...
    'hyperparamsPrior_deg_s02',hyperparamsPrior_deg_s02,...
    'hyperparamsPrecDeg0',hyperparamsPrecDeg0,'hyperparamsPrecSigma0',hyperparamsPrecSigma0,...
    'profileVecSize',profileVecSize,'priorPlotsOn',priorPlotsOn,...
    'regProfilePlotsOn',regProfilePlotsOn,'numRegulatorsInFig',numRegulatorsInFig,...
    'minSwTimeDefault',minSwTimeDefault,'minSwTime',minSwTime) ;

clear i profiles priors nregsPriors regsPriors threshPriors

lastr = num2str(hyperparamsLambda0) ;
lastr(lastr == '.') = ',' ;

if strcmp(thresholdsPriorType,'gradient')
    str1 = ['prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'pct' num2str(thresholdsPriorLBoundStandGrads{2}) 'la' lastr];
    prompt = ['do you wish to use the name: "prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'pct' num2str(thresholdsPriorLBoundStandGrads{2}) 'la' lastr  '" for the output file?(0/1)'] ;
    if input(prompt)
        save([parentpath '/output/prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'pct' num2str(thresholdsPriorLBoundStandGrads{2}) 'la' lastr])
    else
        prompt = 'please type the filename:' ;
        str1=input(prompt,'s');
        save([parentpath '/output/' str1])
    end
elseif strcmp(thresholdsPriorType,'uniform')
    str1 = ['prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'la' lastr] ;
    prompt = ['do you wish to use the name: "prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'la' lastr '" for the output file?(0/1)'] ;
    if input(prompt)
        save([parentpath '/output/prof_priors_prop' namestr '_' targetName '_RegTYPE' regulatorsPriorType '_ThreshTYPE' thresholdsPriorType 'la' lastr])    
    else
        prompt = 'please type the filename:' ;
        str1=input(prompt,'s');
        save([parentpath '/output/' str1])
    end
end
priors_filename = str1;

cd('..'); 
cd('scripts/')

