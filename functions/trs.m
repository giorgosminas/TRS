function output = trs(basenet_filename,priors_filename,targetName,iterations)

currentpath = cd('..');
parentpath = pwd();
cd(currentpath);

load([parentpath '/output/' basenet_filename],'net') ;

regNames = net.netForInterType{strcmp(net.netForInterType{:,3},targetName),1} ;
% regNamesAll = net.regulatorNamesForInterType(strcmp(net.regulatorNamesForInterType(:,1),regNames)) ;

load([parentpath '/output/' priors_filename],'profilesPriors') ;

profilesOutput = profilesPriors.profilesOutput ;
priorsPropOutput = profilesPriors.priorsPropOutput ;

%% mcmc 
proposalsRwSd = 50 ; 
tic; mcmcOutput = createRJMCMC(profilesOutput,priorsPropOutput,'iterations',iterations,'proposalsRwSd',proposalsRwSd) ; telmcmc=toc ;

%% prepare final output
tic ; finalOutput = createFinalOutput(profilesOutput,priorsPropOutput,mcmcOutput,'regulatorNames',regNames) ; telfinalo=toc ;

output.mcmcOutput = mcmcOutput ;
output.finalOutput = finalOutput ;
output.telmcmc = telmcmc ;
output.telfinalo = telfinalo ;

end

