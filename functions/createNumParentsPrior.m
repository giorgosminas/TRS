function numParentsPriorOutput = createNumParentsPrior(type,parentPriorProbs,maxNumParents,numSelParentsPriorLambda)

numParents = length(parentPriorProbs) ;
numParentsPrior = zeros(numParents+1,1) ;

if (strcmp(type,'poisson'))
    
    for z = 0:maxNumParents
        numParentsPrior(z+1) = truncPoisson(z,numSelParentsPriorLambda,maxNumParents) ;
    end
    numParentsPriorOutput.probs = numParentsPrior ;
    numParentsPriorOutput.type = 'poisson' ;
    
elseif (strcmp(type,'rangeProbs'))
    
    numParentsPrior(1) = prod(1 - parentPriorProbs) ;
    
    for z = 1:maxNumParents
        combis = nchoosek(1:numParents,z) ;
        [leCombis,~] = size(combis) ;
        for i=1:leCombis
            combi = combis(i,:) ;
            numParentsPrior(z+1) = numParentsPrior(z+1) + prod(parentPriorProbs(combi))*prod(1 - parentPriorProbs(setdiff(1:numParents,combi))) ;
        end
        
    end
elseif (strcmp(type,'product'))
    numParentsPrior(1) = prod(1 - parentPriorProbs) ;
    
    for z = 1:maxNumParents
        combis = nchoosek(1:numParents,z) ;
        [leCombis,~] = size(combis) ;
        for i=1:leCombis
            combi = combis(i,:) ;
            numParentsPrior(z+1) = numParentsPrior(z+1) + prod(parentPriorProbs(combi))*prod(1 - parentPriorProbs(setdiff(1:numParents,combi))) ;
        end
        
    end
    numParentsPrior = numParentsPrior'.*poisspdf(0:numParents,numSelParentsPriorLambda) ;
    numParentsPrior = numParentsPrior/sum(numParentsPrior) ;
    
end    

numParentsPriorOutput.probs = numParentsPrior ;
numParentsPriorOutput.type = 'rangeProbs' ; 

end



function prob = truncPoisson(k,lambda,kmax)

if (k > -1 && k < kmax+1)
    prob = exp(-lambda)*lambda^k/factorial(k) ;
else
    prob = 0 ;
end

end