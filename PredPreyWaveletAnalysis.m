%% Initialisation
% the predator-prey data analysed here was published by:
%   Blasius, B., Rudolf, L., Weithoff, G., Gaedke, U., and Fussmann, G. F.
%   Long-term cyclic persistence in an experimental predator–prey system.
%   Nature (2020), 577(7789):226–230.
addpath('Tools');
c = '1'; % Number of analyized Experiment
s = ['Data/Blasius2020/C',c,'.csv'];
Data = csvread(s,1,0);
[Data, t] = treatData(Data',1);
prey = Data(2,:);
pred = Data(3,:);
dt= 1;% time steps in days
pad=1;% pad the time series with zeros
nvoices = 100;% total number of voices100
noctave = 3;% total number of octaves
s0=2*dt; % starting scale
mother='MORLET'; %mother wavelet used
param=6; % wavenumber k0 of morlet wavelet
N = 1000; % generate 1000 surrogates 
alpha = 0.05; % significance level
methodSurr = 3; % method to generate surrogates 1= scrambling, 2= AAFT, 3= IAAFT_td then fc= 0.05 (typically)
if methodSurr == 1
    mode = 'Scram';
    fc = 0;
elseif methodSurr == 2
    mode = 'AAFT';
    fc = 0;
elseif methodSurr == 3
    mode = 'IAAFTtd';
    fc = 0.05;
end
%% Estimation of wavelet coherence (waveCo), wavelet spectra of prey (powerX) and predator (powerY) and wavelet cross-spectrum (powerXY)
[waveCo, powerX, powerY, powerXY, period, scale, coi] = waveCohe(prey, pred, dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX3, ~, ~, ~] = waveCohe(prey, Data(4,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX4, ~, ~, ~] = waveCohe(prey, Data(5,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX5, ~, ~, ~] = waveCohe(prey, Data(6,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX6, ~, ~, ~] = waveCohe(prey, Data(7,:), dt, pad, nvoices, noctave, s0, mother, param);
%% Generate surrogate data and determine critical values
m = size(scale,2);
n = size(powerX,2);
surrX = SurrogateData(pred, N, methodSurr,fc);
surrY = SurrogateData(prey, N, methodSurr,fc);
surrX = squeeze(surrX);
surrY = squeeze(surrY);
SurrwaveCo = zeros(N,m,n);
tic
parfor i=1:N
    [SurrwaveCo(i,:,:), ~, ~, ~, ~, ~, ~] = waveCohe(squeeze(surrX(i,:)), squeeze(surrY(i,:)), dt, pad, nvoices, noctave, s0, mother, param);
end
toc
SurrsigCohe = squeeze(quantile(SurrwaveCo,1-alpha));
%% Length of coherent and non-coherent regimes
[LengCoheOsciExp, LengUnCoheExp, dominaF, FreqInd, high, low] = LengthRegimesExp(waveCo,SurrsigCohe,powerXY,nvoices, period, coi);
[LengCoheOsciSurr, LengUnCoheSurr] = LengthRegimesSurr(SurrwaveCo,SurrsigCohe,high,low,period, coi);
%% Plot
figure('Units', 'inches','Position', [0 0 7 (5/7)*7],'renderer','painters')
    set(0,'DefaultTextFontname', 'CMU Serif')
    set(0,'DefaultAxesFontName', 'CMU Serif')
    set(0,'DefaultAxesFontSize', 11.5)
    set(0,'DefaultTextFontSize', 11.5)

%Plot predator prey
s11 = subplot(3,4,4);
plot(prey(50:86)/max(prey),pred(50:86)/max(pred))
text(-0.3,1,'b','Units','normalized','FontSize',13)
xlim([0 1]); ylim([0 1]);
xlabel('Prey'); ylabel('Predator')

% plot a trajecory in phase space
s12 = subplot(3,4,[1 2 3]);
ylim([0 1]);
plot(t,prey/max(prey))
hold on
plot(t,pred/max(pred))
xlim([0 t(length(t))])
text(-0.07,1.02,'a','Units','normalized','FontSize',13)
hold off
ylabel('Abundance');

% plot time mean of wavelet cross spectrum
high = period(high)*ones(1,length(t));
low = period(low)*ones(1,length(t));
s21 = subplot(3,4,8);
plot(log(mean(powerXY,2)), log2(period));
text(-0.15,1.02,'d','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]);
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
xlabel('Power');

 %plot wavelet coherence
s22 = subplot(3,4,[5 6 7]);
pcolor(t,log2(period),waveCo)
colormap parula
axis tight
shading flat
hold on
contour(t,log2(period),waveCo-SurrsigCohe,[0 0],'k','linewidth',1)
plot(t,log2(dominaF),'k','linewidth',2)
plot(t, log2(coi),'k')
plot(t, log2(high),'k')
plot(t, log2(low),'k')
text(-0.07,1.1,'c','Units','normalized','FontSize',13)
hold off
ylim([min(log2(period)) max(log2(period))]);
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
xlabel('Times [d]'); ylabel('Period [d]');

%plot phase relationships
n = size(t,2);
anglePlotXY = zeros(1,n);
anglePlotX3 = zeros(1,n);
anglePlotX4 = zeros(1,n);
anglePlotX5 = zeros(1,n);
timp = zeros(1,n);
tunimp = zeros(1,n);
for i = 1:n
    if  FreqInd(i) > 1
        anglePlotXY(i) = angle(powerXY(FreqInd(i),i));
        anglePlotX3(i) = angle(powerX3(FreqInd(i),i));
        anglePlotX4(i) = angle(powerX4(FreqInd(i),i));
        anglePlotX5(i) = angle(powerX5(FreqInd(i),i));
    end
end
anglePlotXYi = zeros(1,n);
anglePlotX3i = zeros(1,n);
anglePlotX4i = zeros(1,n);
anglePlotX5i = zeros(1,n);
k=0;
for i = 1:n
    if  anglePlotXY(i) ~= 0
        k = k +1;
        anglePlotXYi(k) = anglePlotXY(i);
        anglePlotX3i(k) = anglePlotX3(i);
        anglePlotX4i(k) = anglePlotX4(i);
        anglePlotX5i(k) = anglePlotX5(i);
    end
end
anglePlotXYi = anglePlotXYi(1:k);
anglePlotX3i = anglePlotX3i(1:k);
anglePlotX4i = anglePlotX4i(1:k);
anglePlotX5i = anglePlotX5i(1:k);
s61del =  subplot(3,4,12);
s31 = polaraxes('Units',s61del.Units,'Position',s61del.Position);
delete(s61del);
polarhistogram(anglePlotXYi,9,'Normalization','pdf')
hold on
polarhistogram(anglePlotX3i,9,'Normalization','pdf')
polarhistogram(anglePlotX4i,9,'Normalization','pdf')
polarhistogram(anglePlotX5i,9,'Normalization','pdf')
text(-0.17,1.1,'f','Units','normalized','FontSize',13)
hold off

% plot critical of the wavelet coherence
s32 = subplot(3,4,[9 10 11]);
pcolor(t,log2(period),SurrsigCohe);
text(-0.07,1.02,'e','Units','normalized','FontSize',13)
colormap parula
axis tight
shading flat
hold on
plot(t, log2(coi),'k')
ylim([min(log2(period)) max(log2(period))]);
hold off
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
ylabel('Period [d]'); xlabel('Time [d]'); 
    
s11.Position = s11.Position + [0.03 0.04 -0.03 -0.07];
s12.Position = s12.Position + [0 0 0.01 -0.03];
s21.Position = s21.Position + [0 0.05 0 -0.03];
s22.Position = s22.Position + [0 0.05 0.01 -0.03];
s31.Position = s31.Position + [0 0 -0.01 -0.01];
s32.Position = s32.Position + [0 0.05 0.01 -0.03];
s = ['Figures/ExpC',c,mode];        
export_fig(s, '-pdf', '-transparent');
