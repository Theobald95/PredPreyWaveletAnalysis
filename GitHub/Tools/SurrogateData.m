function surrData = SurrogateData(x, N, method, fc)
% x: vector of time series
% N: Number of surrogates
% method: Randomization method to generate surrogates 1=scrambling, 2=AAFT,
% 3=IAAFT_td
% fc: cutoff frequency for IAAFT_td method (typical choice 0.05)

if method == 1
    surrData = ScramData(x, N);
elseif method == 2
    surrData = AAFTData(x, N);
elseif method == 3
    surrData = IAAFTtdData(x, N, 0.5-fc);
end
end
%%
function surrData = ScramData(x, N)
n = size(x,2);
surrData = zeros(N,n);
for i=1:N
    
    surrData(i,:) = x(1,randperm(n));
end
end
%%
function surrData = AAFTData(x, N)
y = gaussData(x);
n1 = size(x,2);
y = TaperData(y);
y = PaddData(y);
n = size(y,2);
m = size(y,1);
surrData= zeros(N,m,n);
parfor i = 1:m
    X = fft(y(i,:));
%     rng('shuffle')
    randPhases = rand(N,n/2-1)*2*pi;
    randPhasesSym = [zeros(N,1), randPhases, zeros(N,1),-fliplr(randPhases)];
    surrData(:,i,:) = X.*exp(1i*randPhasesSym);
    surrData(:,i,:) = ifft(squeeze(surrData(:,i,:)),[],2,'symmetric');
end
surrData = surrData(:,:,1:n1);
parfor i= 1:m
    for j = 1:N
        surrData(j,i,:) = invgaussData(x,squeeze(surrData(j,i,:)));
    end
end
end
%%
function surrData = IAAFTtdData(x, N, fc)
n1 = length(x);
x = PaddData(x);
n = length(x);

fc = floor(n*fc);
fr = n/2-fc;
t = 0:n-1;
stdval=std(x);
meanx = mean(x);
[xsorted,~] = sort(x-meanx);
g  = (((x(end)-meanx)-(x(1)-meanx))/(n-1))*(t-1);
xdetrend = x-meanx -g;
surrData = zeros(N,n1);
X = fft(xdetrend);
sdtmodules = abs(X);
sdtphase = angle(X);
for i = 1:N
%     rng('shuffle')
    randPhases = rand(1,n/2)*2*pi;
    randPhasesSym = [ randPhases, -fliplr(randPhases)];
    xdetrendRP = X.*exp(1i*randPhasesSym);
    xdetrendRP = ifft(squeeze(xdetrendRP),[],2,'symmetric');
    xtrendRP =xdetrendRP + g;
    xtrendprimeprime  = xtrendRP;
    accerror=0.0001;
    amperror(1)=100;
    specerror(1)=100;   
    counter=1;
    while (amperror(counter) > accerror) && (specerror(counter) > accerror)
        xtrendprime = xtrendprimeprime;
        [~,shuffind]=sort(real(xtrendprime));
        xtrendprimeprime(shuffind)=xsorted;
        ampdiff=mean(abs(real(xtrendprimeprime)-real(xtrendprime)));
        amperror(counter+1) = ampdiff/stdval;
        
        xtrendprime = xtrendprimeprime;
        xdetrendprime = xtrendprime - g;
        Xprime = fft(xdetrendprime);
        sdtphaseprime = angle(Xprime);
        phase = [sdtphase(1:fr), sdtphaseprime(fr+1:fr+2*fc), sdtphase(fr+2*fc+1:n)];
        xdetrendprimeprime = ifft(sdtmodules.*exp(1i*phase));
        xtrendprimeprime = xdetrendprimeprime +g;
        specdiff=mean(mean(abs(real(xtrendprimeprime))-real(xtrendprime)));
        specerror(counter+1) = specdiff/stdval;
        
        toterror=amperror(counter+1)+specerror(counter+1);
        oldtoterr=amperror(counter)+specerror(counter);
        if (oldtoterr-toterror)/toterror < (accerror/10)
            amperror(counter+1)=-1;
        end
        
        counter=counter+1;
    end
    surrData(i,:) = xtrendprimeprime(1:n1)+meanx;
end
end
%%
function X = gaussData(x)
m = size(x,1);
n = size(x,2);
X = randn(m,n);
for i= 1:m
    [~, indexx] = sort(x(i,:),2);
    xmu = mean(x(i,:));
    xstd = std(x(i,:));
    sortedX = sort(X(i,:),2);
    X(i,indexx) = sortedX*xstd + xmu;
end
end
%%
function X = invgaussData(x,s)
[~, indexs] = sort(s);
X = x;
sortedX = sort(x);
X(indexs) = sortedX;
end
%%
function paddData = PaddData(X)
N = length(X);
l=log2(N);
p=ceil(l);  
N2 = 2^p;
paddData = zeros(size(X,1),N2);
paddData(:,1:N)= X;
end
%%
function tapData = TaperData(X)
N = size(X,2);
barlett = 1:N+1;
barlett = 1-abs((barlett-0.5*(N+1))/(0.5*(N+1)));
tapData = zeros(size(X));
for i = 1:size(X,1)
    tapData(i,:)= X(i,:).*barlett(1:N);
end
tapData = (N/ sum(barlett(1:N)))*tapData;
end