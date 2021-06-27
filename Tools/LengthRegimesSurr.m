function [LengCoheOsciSurr, LengUnCoheSurr] = LengthRegimesSurr(SurrwaveCo,SurrsigCohe,high,low, period, coi)
% SurrwaveCo: wavelet coherence of surrogate data
% SurrsigCohe: critical value for wavelet coherence
% high: highest period to consider, low: lowest period to consider
% period: period-length of oscilattions
% coi: cone of influence
% returns the coherent oscillations regime and non-coherent oscillation regime length for the experiment
N = size(SurrwaveCo,1);
n = size(SurrwaveCo,3);
LengCoheOsciSurr = cell(N,1);
LengUnCoheSurr = cell(N,1);
for i=1:N
    LengCoheOsciSurr{i} = [];
    LengUnCoheSurr{i} = [];
    FreqIndShuff = DomFreqShuff(squeeze(SurrwaveCo(i,:,:)),SurrsigCohe,high,low);
    r = 0;
    l = 0;
    for j = 1:n
        if FreqIndShuff(j) == high -1 || FreqIndShuff(j) == low
            FreqIndShuff(j) = 1;
        end
        if period(FreqIndShuff(j))> coi(j)
            FreqIndShuff(j) = 1;
        end
        if FreqIndShuff(j) > 1
            if j > 2 && FreqIndShuff(j-1) > 1
                LengCoheOsciSurr{i}(r) = LengCoheOsciSurr{i}(r) +1;
            else
                r = r+1;
                LengCoheOsciSurr{i}(r)= 1;
            end
        else
            if j> 2 && FreqIndShuff(j-1) == 1
                LengUnCoheSurr{i}(l) = LengUnCoheSurr{i}(l) +1;
            else
                l = l+1;
                LengUnCoheSurr{i}(l)= 1;
            end
        end
    end
%             temp = [LengCoheOsciSurr{i}(:)];
%             temp = sum(temp'.^2)/sum(temp');
%             TShuff(i)=temp;
end
LengCoheOsciSurr = [LengCoheOsciSurr{:}];
LengUnCoheSurr = [LengUnCoheSurr{:}];
end
%%
function FreqInd = DomFreqShuff(Cohe,sigCohe,high,low)
n = size(Cohe,2);
FreqInd = ones(1,n);
for i = 1:n
    [~, index] = max(Cohe(low:high,i));
    if Cohe(low+index-1,i) > sigCohe(low+index-1,i)
        FreqInd(i) = low+index-1;
    end
end
end
