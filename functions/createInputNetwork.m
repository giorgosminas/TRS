function net = createInputNetwork(net_filename,interactionType,names_filename) 

% creates the network mat files based on the csv files of the network
% the network mat files contain a structure that gives the
% names of regulators and target in the csv file
% and the interactions of type interactionType 

currentpath = cd('..');
parentpath = pwd();
cd(currentpath);

netAll = readtable([parentpath '/input/' net_filename]) ;
namesAll = readtable([parentpath '/input/' names_filename]) ;
ninteractions = height(netAll) ;

interactionTypeTrue = zeros(ninteractions,1) ;
for i=1:ninteractions
    if strcmp(netAll{i,2},interactionType)
        interactionTypeTrue(i) = 1 ;
    end
end

netForInterType = netAll(interactionTypeTrue==1,:) ;

targetNamesForInterType = unique(netForInterType(:,3),'stable') ;
regulatorNamesForInterType = unique(netForInterType(:,1),'stable') ;

targetNamesForInterType1 = cell(height(targetNamesForInterType),width(namesAll)) ;
for i=1:height(targetNamesForInterType)
    idx = strcmp(namesAll{:,1},targetNamesForInterType{i,1}) ;
    for j=1:width(namesAll)
        targetNamesForInterType1(i,j) = namesAll{idx,j} ;
    end
end
targetNamesForInterType = targetNamesForInterType1 ;
clear targetNamesForInterType1 

regulatorNamesForInterType1 = cell(height(regulatorNamesForInterType),width(namesAll)) ;
for i=1:height(regulatorNamesForInterType)
    idx = strcmp(namesAll{:,1},regulatorNamesForInterType{i,1}) ;
    for j=1:width(namesAll)
        regulatorNamesForInterType1(i,j) = namesAll{idx,j} ;
    end
end
regulatorNamesForInterType = regulatorNamesForInterType1 ;
clear regulatorNamesForInterType1 

net = struct ;
net.allnames = namesAll;
net.allinteractions = netAll;
net.netForInterType = netForInterType ;
net.regulatorNamesForInterType =regulatorNamesForInterType;
net.targetNamesForInterType = targetNamesForInterType;


end