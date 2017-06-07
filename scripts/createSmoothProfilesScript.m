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

% clear all  

% run from directory '.../trRegCode/Scripts'
cd('..'); 
parentpath = pwd();
cd('functions/')

global dataset_filenames net_filename smoothProfiles_filename namestr
if isempty(dataset_filenames) || isempty(net_filename)
    dlg_title = 'Filenames & Delay time' ;
    prompt = {'data filename','data filename','data filename','data filename','data filename','data filename',...
        'Base network filename','Delay time'} ;
    num_lines = 1 ;
    defaultPreC_dlg = {'','','','','','',''} ;
    precomputePriorsdlg = inputdlg(prompt,dlg_title,num_lines,defaultPreC_dlg) ;
    dataset_filenames=[];
    for i=1:length(prompt)-1
        if ~isempty(precomputePriorsdlg{i})
            dataset_filenames{i}=precomputePriorsdlg{i};
        end
    end
    net_filename = precomputePriorsdlg{end-1};
    delayTime = precomputePriorsdlg{end} ;
end

delayTime = 1 ;
nBootSamples = 100 ;
smoothPar = [0.1 0.1] ; % smoothing parameters. The first is for the first smoothing spline and the second for the bootstrap samples of smoothing splines.

smoothProfiles = createSmoothProfiles(net_filename,dataset_filenames,'delayTime',delayTime,'nBootSamples',nBootSamples,'smoothPar',smoothPar) ;

str = [] ;
for i=1:length(dataset_filenames)
    str = [str dataset_filenames{i}(5:end-4)] ;
end
namestr = str ; 

clear i        
save([parentpath '/output/smoothProfiles' str])
smoothProfiles_filename = ['smoothProfiles' str];

cd('..'); 
cd('scripts/')
