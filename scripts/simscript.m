cd1 = cd ;
mypathtoTRS = '/Users/giorgosminas/Google Drive/Topics/PRESTA/Network modelling/Code/TRSNEW1' ;
cd([mypathtoTRS '/scripts'])
%%
% If running this script multiple times, run the next lines once to create dataset files and then "comment" them
global infodata_filename dataset_filename dataset_filenames
infodata_filename = 'info_sim1.mat' ;
dataset_filename = 'data_sim1.csv' ;
createDatasetScript
infodata_filename = 'info_sim2.mat' ;
dataset_filename = 'data_sim2.csv' ;
createDatasetScript
dataset_filenames = {'data_sim1.mat', 'data_sim2.mat'} ;

%%
% Similarly, run the next lines once to create the network files and then "comment" them
global net_filename interactionType names_filename
net_filename = 'basenetSim.csv' ;
interactionType = 'hyp' ;
names_filename = 'namesSim.csv' ;
createNetworkScript

%%
% Similarly, run the next lines to create the smooth profiles file and then "comment" them
createSmoothProfilesScript

%%
% Similarly, run the next lines to create the priors and proposals file and then "comment" them unless you want to create different priors
global targetName 
targetName = 'T';
createProfilesPriorsAndProposalDistrsScript

%%
% the next line calls the main script. Keep it "uncommented"
trsScript

cd(cd1)

