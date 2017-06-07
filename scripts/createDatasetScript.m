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

global infodata_filename dataset_filename

if isempty(infodata_filename)
    prompt='please provide the name of the information mat file (format: filename.mat, see examples in input folder)';
    infodata_filename = input(prompt) ; % mat file with various variables with information about the data
end
if isempty(dataset_filename)
    prompt='please provide the name of the dataset csv file (format: filename.csv, see examples in input folder)';
    dataset_filename = input(prompt) ; % csv file with data of one experiment
end

dataset = createDataset(infodata_filename,dataset_filename) ;

dataset_filename = [dataset_filename(1:end-4) '.mat'] ;
save([parentpath '/output/' dataset_filename(1:end-4)])

cd('..'); 
cd('scripts/')

%%% RUN THIS SCRIPT FOR ALL DATASETS, altering the 
% infodata_filename and dataset_filename
