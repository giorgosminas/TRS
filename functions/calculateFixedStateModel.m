function [initialExps,birthRates,XAll] = calculateFixedStateModel(timepoints,swTimes,states,statesUsed,degRate,YAll,ntimepoints,replicates,lsquares)

global options 
regFunc = @fixedStateSwitchModelExact;

XAll = zeros(sum(ntimepoints.*replicates),length(swTimes) + length(statesUsed));

currX = 1;
for m = 1:length(swTimes)
    X = regFunc(timepoints{m},swTimes{m},states{m},degRate,statesUsed);
    for z = 1:replicates(m)
        XAll(currX:(currX-1+ntimepoints(m)),m) = X(:,1);
        XAll(currX:(currX-1+ntimepoints(m)),(length(swTimes)+1):end) = X(:,2:end);
        
        currX = currX + ntimepoints(m);
    end
end

if (strcmp(lsquares.type,'wls'))
    Qinvroot = lsquares.QAll^(-1/2) ;
    XAll = Qinvroot*XAll ; 
    YAll = Qinvroot*YAll ;
end
alpha = (XAll' * XAll) \ XAll' * YAll ;
% noFail = sum(alpha < 0) == 0 ;
% if ~noFail
%     alpha = lsqlin(XAll,YAll,[],[],[],[],zeros(size(XAll,2),1),[],alpha,options) ;  
% end

initialExps = alpha(1:length(swTimes));
birthRates = alpha((length(swTimes)+1):end);

end


function [X] = fixedStateSwitchModelExact(timescale,swTimes,states,deg_m,statesUsed)

n = length(timescale) ;
X = zeros(n,length(statesUsed)+1);

X(:,1) = exp(-deg_m.*timescale) ;  % the first column corresponds to M(O)

taub = zeros(length(timescale),1);
taub(1) = states(1);
for p = 2:length(timescale)    
    idx = find(swTimes > timescale(p),1);
    if isempty(idx)
        taub(p) = states(end);
    else
        taub(p) = states(idx);
    end
end

[timeresolution,sIdx] = sort([timescale swTimes]);
taub = [taub' states(2:end)];
taub = taub(sIdx);

%faster to use diff
db = diff(timeresolution);
timeresolution = [timeresolution(db > 0) timeresolution(end)];
taub = [taub(db > 0) taub(end)];

xx1 = exp(deg_m*timeresolution);
xx = (1/deg_m)*(xx1(2:end) - xx1(1:end-1));

%trying to vectorise the code for speed
et = exp(-deg_m .* timescale)'; %replaces eti call

[~,timeIdx] = ismember(timescale,timeresolution); %replaces the find(timeresolution=ti) call

for k=1:length(statesUsed)
    s = statesUsed(k);
    for i=2:n
        xik = sum(xx(taub(1:(timeIdx(i)-1)) == s));
        X(i,k+1) = xik;
    end
    X(:,k+1) = X(:,k+1) .* et;
end

end

