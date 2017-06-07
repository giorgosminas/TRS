function output = createRegulatorsPrior(rangeReg)  

nregulators = size(rangeReg,2) ;

mur = mean(rangeReg) ; % mean range of each reg across Bootstrap samples
stdr = std(rangeReg) ; % std of range of each reg across Bootstrap samples

coeffr = mur./stdr ; % inverse of coefficient of variation of range of each reg across Bootstrap samples

regPriorsMatrixAllcomb = cell(1,min(8,nregulators)) ;

p1 = normcdf(coeffr,mean(coeffr),std(coeffr))' ; %*sigmaFactor
p1 = p1/sum(p1) ;

regPriorsMatrixAllcomb{1} = p1 ;
clear rangeReg

if (nregulators >= 2)
    
    combiCn = nchoosek(nregulators,2) ;
    combin = nchoosek(1:nregulators,2) ;
    regPriorsMatrixAllcomb{2} = zeros(combiCn,2+1) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{2} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(2)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    % this gets the priors of each combination
    % each dimension of selected parents gets into a cell
    regPriorsMatrixAllcomb{2} = regPriorsMatrixAllcombn ;

end

if (nregulators >= 3)
    % the third is a 3d matrix with (j1,j2,i1) entry the prob of pick i1 when j1,
    % j2 are already picked
    % note that (j1,i1,i1), (i1,j2,i1) and (i2,i2,i1) are all zero
    
    combiCn = nchoosek(nregulators,3) ;
    combin = nchoosek(1:nregulators,3) ;
    regPriorsMatrixAllcomb{3} = zeros(combiCn,3+1) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{3} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(3)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1-p1(pe(q,1))-p1(pe(q,2)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    regPriorsMatrixAllcomb{3} = regPriorsMatrixAllcombn ;
end

if (nregulators >= 4)
    % the fourth is a 4d matrix with (j1,j2,j3,i1) entry the prob of pick i1 when j1,
    % j2,j3 are already picked
    % note that (j1,j2,i1,i1), (j1,i1,j3,i1), (i1,j2,j3,i1), (i2,i2,j3,i1), (i2,j2,i2,i1), (j1,i2,i2,i1) are all zero
    
    combiCn = nchoosek(nregulators,4) ;
    combin = nchoosek(1:nregulators,4) ;
    regPriorsMatrixAllcomb{4} = zeros(combiCn,5) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{4} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(4)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1- p1(pe(q,1)) - p1(pe(q,2)))*p1(pe(q,4))/(1-p1(pe(q,1)) - p1(pe(q,2)) - p1(pe(q,3)))   ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:5) = combin(i,:) ;
    end    
    
    regPriorsMatrixAllcomb{4} = regPriorsMatrixAllcombn ;
end

if (nregulators >= 5)
    % the 5th is a 5d matrix with (j1,j2,j3,j4,i) entry being the pr
    % of picking i when j1,j2,j3,j4 are picked
    % this is p1(i)/(1-p1(j1)-p1(j2)-p1(j3)-p1(j4))
    % note that lots of entries will be 0 to avoid double selection
    % we first need to create a 5d matrix with all combination of sums
    
    
    combiCn = nchoosek(nregulators,5) ;
    combin = nchoosek(1:nregulators,5) ;
    regPriorsMatrixAllcomb{5} = zeros(combiCn,1+5) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{5} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(5)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1-p1(pe(q,1)) - p1(pe(q,2)))*p1(pe(q,4))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3)))*p1(pe(q,5))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    regPriorsMatrixAllcomb{5} = regPriorsMatrixAllcombn ;
end


if (nregulators >= 6)
    % the 6th is a 6d matrix with (j1,j2,j3,j4,j5,i) entry being the pr
    % of picking i when j1,j2,j3,j4,j5 are picked
    % this is p1(i)/(1-p1(j1)-p1(j2)-p1(j3)-p1(j4))
    % note that lots of entries will be 0 to avoid double selection
    % we first need to create a 4d matrix with all combination of sums
        
    combiCn = nchoosek(nregulators,6) ;
    combin = nchoosek(1:nregulators,6) ;
    regPriorsMatrixAllcomb{6} = zeros(combiCn,1+6) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{6} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(6)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1-p1(pe(q,1)) - p1(pe(q,2)))*p1(pe(q,4))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3)))*p1(pe(q,5))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4)))*p1(pe(q,6))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    regPriorsMatrixAllcomb{6} = regPriorsMatrixAllcombn ;

end

if (nregulators >= 7)
    % the 6th is a 6d matrix with (j1,j2,j3,j4,j5,i) entry being the pr
    % of picking i when j1,j2,j3,j4,j5 are picked
    % this is p1(i)/(1-p1(j1)-p1(j2)-p1(j3)-p1(j4))
    % note that lots of entries will be 0 to avoid double selection
    % we first need to create a 4d matrix with all combination of sums
    
    combiCn = nchoosek(nregulators,7) ;
    combin = nchoosek(1:nregulators,7) ;
    regPriorsMatrixAllcomb{7} = zeros(combiCn,1+7) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{7} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(7)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1-p1(pe(q,1)) - p1(pe(q,2)))*p1(pe(q,4))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3)))*p1(pe(q,5))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4)))*p1(pe(q,6))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5)))*p1(pe(q,7))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5))-p1(pe(q,6)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    regPriorsMatrixAllcomb{7} = regPriorsMatrixAllcombn ;

end

if (nregulators >= 8)
    % the 6th is a 6d matrix with (j1,j2,j3,j4,j5,i) entry being the pr
    % of picking i when j1,j2,j3,j4,j5 are picked
    % this is p1(i)/(1-p1(j1)-p1(j2)-p1(j3)-p1(j4))
    % note that lots of entries will be 0 to avoid double selection
    % we first need to create a 4d matrix with all combination of sums
   
    combiCn = nchoosek(nregulators,8) ;
    combin = nchoosek(1:nregulators,8) ;
    regPriorsMatrixAllcomb{8} = zeros(combiCn,1+8) ;
    regPriorsMatrixAllcombn = regPriorsMatrixAllcomb{8} ;
    for i=1:combiCn
        p = 0 ;
        pe = perms(combin(i,:)) ;
        for q = 1:factorial(8)
            p = p + p1(pe(q,1))*p1(pe(q,2))/(1-p1(pe(q,1)))*p1(pe(q,3))/(1-p1(pe(q,1)) - p1(pe(q,2)))*p1(pe(q,4))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3)))*p1(pe(q,5))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4)))*p1(pe(q,6))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5)))*p1(pe(q,7))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5))-p1(pe(q,6)))*p1(pe(q,8))/(1-p1(pe(q,1))-p1(pe(q,2))-p1(pe(q,3))-p1(pe(q,4))-p1(pe(q,5))-p1(pe(q,6))-p1(pe(q,7)))  ;
        end
        regPriorsMatrixAllcombn(i,1) = p ;
        regPriorsMatrixAllcombn(i,2:end) = combin(i,:) ;
    end
    
    regPriorsMatrixAllcomb{8} = regPriorsMatrixAllcombn ;

end

output.regPriorsMatrixAllcomb = regPriorsMatrixAllcomb ;

end





