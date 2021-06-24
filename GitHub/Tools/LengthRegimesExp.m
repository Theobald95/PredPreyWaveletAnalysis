function [LengCoheOsciExp, LengUnCoheExp, dominaF, FreqInd, high, low] = LengthRegimesExp(waveCo,sigCohe,powerXY,nvoices, period, coi)
n = size(powerXY,2);
LengCoheOsciExp =  [];
LengUnCoheExp = [];
[dominaF, FreqInd, high, low] = DomFreq(waveCo,sigCohe,powerXY,nvoices, period);
r = 0;
l = 0;
for j = 1:n
    if FreqInd(j) == high -1 || FreqInd(j) == low
            FreqInd(j) = 1;
            dominaF(j) = 0;
    end
    if period(FreqInd(j))> coi(j)
        FreqInd(j) = 1;
        dominaF(j) = 0;
    end
    if FreqInd(j) > 1
        if j > 2 && FreqInd(j-1) > 1
            LengCoheOsciExp(r) = LengCoheOsciExp(r) +1;
        else
            r = r+1;
            LengCoheOsciExp(r)= 1;
        end
    else
        if j> 2 && FreqInd(j-1) == 1
            LengUnCoheExp(l) = LengUnCoheExp(l) +1;
        else
            l = l+1;
            LengUnCoheExp(l)= 1;
        end
    end
end
end
%%
function [dominaF, FreqInd,high,low] = DomFreq(Cohe,sigCohe,powerX, nvoices, scale)
n = size(Cohe,2);
[~, mid] = max(mean(powerX,2));
plumin = round(0.6*nvoices);
high = mid + plumin;
low = mid - plumin;
if high > size(Cohe,1)
    high = size(Cohe,1);
end
if low < 1
    low = 1;
end
dominaF = zeros(1,n);
FreqInd = ones(1,n);
for i = 1:n
    [~, index] = max(Cohe(low:high,i));
    if Cohe(low+index-1,i) > sigCohe(low+index-1,i)
        dominaF(1,i) = scale(low+index-1);
        FreqInd(i) = low+index-1;
    end
end
end