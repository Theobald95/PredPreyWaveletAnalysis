function [xNeu, tNeu] = treatData(x,dt)
n = size(x,2);
m = size(x,1);
tend = round(x(1,n));
tNeu = 0:dt:tend;
xNeu = zeros(m,length(tNeu));
for k = 2:m
    j= 0;
    xNan = zeros(1,n);
    tNan = zeros(1,n);
    for i = 1:n
            if ~isnan(x(k,i))
                j = j +1;
                xNan(j) = x(k,i);
                tNan(j) = x(1,i);
            end
    end
    xNan = xNan(1:j);
    tNan = tNan(1:j);
    xNeu(k,:) = spline(tNan,xNan(:),tNeu);
    if min(xNeu(k,:)) < 0 
        xNeu(k,:) = xNeu(k,:) - min(xNeu(k,:));
    end
end
end