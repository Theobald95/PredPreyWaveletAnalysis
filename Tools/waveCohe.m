%%
function [waveCo, powerX, powerY, powerXY, period, scale, coi] = waveCohe(x, y, dt, pad, nvoices, noctave, s0, mother, param)
% x: vector of time series 1
% y: vector of time series 2
% dt: sampling interval, pad: 0=don't padd data; 1= padd&taper data
% noctave: number of powers of two, nvoices: number of steps inbetween
% successive powers of 2
% s0: starting scale (typically s0=2*dt), mother: mother wavelet (here:
% mother='MORLET';), param= for 'MORLET' is k0(wavenumber) default 6
        
[waveX, ~, ~, ~] = cwt(x,dt,pad, nvoices, noctave, s0, mother, param);
[waveY, period, scale, coi] = cwt(y,dt,pad, nvoices, noctave, s0, mother, param);
powerX = smoothWave(waveX.*conj(waveX),scale,nvoices);
powerY = smoothWave(waveY.*conj(waveY),scale,nvoices);
powerXY = smoothWave(waveX.*conj(waveY),scale,nvoices);
waveCo = abs(powerXY)./(sqrt(powerX).*sqrt(powerY));
powerX = invTaper(powerX);
powerY = invTaper(powerY);
powerXY = invTaper(powerXY);
end
%%
function [wave, period, scale, coi] = cwt(x,dt,pad, nvoices, noctave, s0, mother, param)
    dj = 1 / nvoices;
    j1 = noctave * nvoices;
    x = x - mean(x);  
    x = x ./ std(x);
    [wave, period, scale, coi] = wavelet(x,dt,pad,dj,s0,j1,mother,param);
end
%%
function smWave = smoothWave(power,scale,nvoices)
[nscale, n] = size(power);
smWave = zeros(nscale,n);
sw_tot = round(0.6*nvoices);
parfor i = 1:nscale
    smWave(i,:) = smoothdata(power(i,:),'gaussian',5*scale(i)); %std of gaussian is 1/5 of window size and you can't change that
end
parfor i = 1:n
    smWave(:,i) = smooth(smWave(:,i), sw_tot);
end
end
%%
function invTap = invTaper(X)
N = size(X,2);
barlett = 1:N+1;
barlett = 1-abs((barlett-0.5*(N+1))/(0.5*(N+1)));
invTap = zeros(size(X));
for i = 1:size(X,1)
    invTap(i,:)= X(i,:)./barlett(1:N);
end
end