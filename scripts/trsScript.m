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

rng('shuffle') 

tic; 
% run from directory '.../trRegCode/Scripts'
cd('..'); 
parentpath = pwd();
cd('functions/')

global dataset_filenames net_filename targetName priors_filename
if isempty(dataset_filenames) || isempty(net_filename) || isempty(targetName) || isempty(priors_filename)
    dlg_title = 'Filenames and Target' ;
    prompt = {'data filename','data filename','data filename','data filename','data filename','data filename',...
        'Base network filename','priors_filename','targetName'} ;
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
    priors_filename = precomputePriorsdlg{end-1};
    targetName = precomputePriorsdlg{end};
end
iterations = 1000;

output = trs(net_filename,priors_filename,targetName,iterations) ;

str = priors_filename(17:end);     
tEl=toc;
save([parentpath '/output/finalOutput' str])


cd('..'); 
cd('scripts/')
