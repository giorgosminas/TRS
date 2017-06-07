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

% run from directory '.../trRegCode/Scripts'
cd('..'); 
parentpath = pwd();
cd('functions/')

global net_filename interactionType names_filename
if isempty(net_filename) || isempty(names_filename)
    prompt='please provide the name of the base network csv file (format: filename.csv, see examples in input folder)';
    net_filename = input(prompt) ; %file with base network structure. It may contain multiple targets & interaction types
end
if isempty(interactionType)
    prompt='please provide the interaction type. This must be one of those interactions in the network file' ;
    interactionType = input(prompt) ;  % type of base interactions
end
if isempty(names_filename)
    prompt='please provide the name of the file containing the names of the regulators and targets in the network (format: filename.csv, see examples in input folder)' ;
    names_filename = input(prompt) ;  % file with (different) names of each regulator & target
end

net = createInputNetwork(net_filename,interactionType,names_filename); 

net_filename = [net_filename(1:end-4) interactionType] ;
save([parentpath '/output/' net_filename])

cd('..'); 
cd('scripts/')
