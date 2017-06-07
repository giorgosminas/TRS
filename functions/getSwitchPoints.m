function [sw,st,activation,switches,states] = getSwitchPoints(activation,timescale,deltaSw)

% keyboard
swPos = find(abs(activation(2:end) - activation(1:(end-1))))+1;
switches = timescale(swPos);
states = [activation(1) activation(swPos)'];

%remove any switch/state where the interval is smaller than the given time
intervals = [switches timescale(end)] - [timescale(1) switches];

sw = switches ; % fliplr(switches);
st = states ; %fliplr(states);
i = intervals ; % fliplr(intervals);
idx = 1:1:length(switches) ; % length(switches):-1:1;
[v,x] = min(i); % 

% if v< deltaSw
%     keyboard
% end

while v < deltaSw
%     
    if x == 1        
        sw = sw(2:end);
        st = st(2:end) ; % [st(1) st(3:end)];
        idx = idx(2:end);
          
    elseif x == length(i) % length(i)
        sw = sw(1:(end-1));
        st = st(1:(end-1)); % [st(1:(end-2)) st(end)];
        idx = idx(1:(end-1));
%         activation(swPos(end):end) = activation(swPos(end)-1) ;
    else
        sw = [sw(1:(x-1)) sw((x+1):end)];
        st = [st(1:(x-1)) st((x+1):end)];
        idx = [idx(1:(x-1)) idx((x+1):end)];
%         activation(swPos(x-1):swPos(x)-1) = activation(swPos(x-1)-1) ;
    end
    i = [sw timescale(end)] - [timescale(1) sw];
% [timescale(end) sw] - [sw timescale(1)];
    
    [v,x] = min(i);
    
end

stf = st ;
swf = sw ;
idx1=1 ;
while idx1<=length(stf)-1
    if stf(idx1+1)-stf(idx1)==0
        stf = stf([1:idx1-1 idx1+1:end]) ;
        swf = swf([1:idx1-1 idx1+1:end]) ;
    else
        idx1 = idx1+1;
    end
end

sw = swf ;
st = stf ;

% sw = fliplr(sw);
% st = fliplr(st);
% i = fliplr(i);
% idx = fliplr(idx);

end

