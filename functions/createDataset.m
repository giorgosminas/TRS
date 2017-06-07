function dataset = createDataset(infodata_filenames,data_filenames)

% creates mat file with structure containing the data and related information

currentpath = cd('..');
parentpath = pwd();
cd(currentpath);

load([parentpath '/input/' infodata_filenames],'info') ;
dataset = info ;
dataset.data = importdata([parentpath '/input/' data_filenames]) ;


end