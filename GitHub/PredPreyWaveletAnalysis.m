%% Initialisation
addpath('Tools');
c = '1'; % Number of analyized Experiment
s = ['Data/Blasius2020/C',c,'.csv'];
Data = csvread(s,1,0);
[Data, t] = treatData(Data',1);
prey = Data(2,:);
pred = Data(3,:);
dt= 1;% time steps in days
pad=1;% pad the time series with zeros
nvoices = 1000;% total number of voices100
noctave = 3;% total number of octaves
s0=2*dt; % starting scale
mother='MORLET'; %mother wavelet used
param=6; % wavenumber k0 of morlet wavelet
N = 100; % generate 1000 surrogates 
alpha = 0.05; % significance level
methodSurr = 2; % method to generate surrogates 1= scrambling, 2= AAFT, 3= IAAFT_td then fc= 0.05 (typically)
if methodSurr == 1
    mode = 'Scram';
elseif methodSurr == 2
    mode = 'Surr';
elseif methodSurr == 3
    mode = 'IAAFTtd';
end
%% Estimation of wavelet coherence and wavelet cross spectra
[waveCo, powerX, powerY, powerXY, period, scale, coi] = waveCohe(prey, pred, dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX3, ~, ~, ~] = waveCohe(prey, Data(4,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX4, ~, ~, ~] = waveCohe(prey, Data(5,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX5, ~, ~, ~] = waveCohe(prey, Data(6,:), dt, pad, nvoices, noctave, s0, mother, param);
[~, ~, ~, powerX6, ~, ~, ~] = waveCohe(prey, Data(7,:), dt, pad, nvoices, noctave, s0, mother, param);
%% Generate surrogate data and determine significance levels
m = size(scale,2);
n = size(powerX,2);
surrX = SurrogateData(pred, N, methodSurr);
surrY = SurrogateData(prey, N, methodSurr);
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
 figure('Units', 'inches','Position', [0 0 6.2 (6.5/5)*6.20],'renderer','painters')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
set(0,'DefaultAxesFontSize', 11.5)
set(0,'DefaultTextFontSize', 11.5)

%Plot predator prey time series
s11= subplot(6,4,4);
plot(prey(50:86)/max(prey),pred(50:86)/max(pred))
text(-0.45,1,'b','Units','normalized','FontSize',13)
xlim([0 1]); ylim([0 1]);
xlabel('Prey','Interpreter','latex'); ylabel('Predator','Interpreter','latex')
s12= subplot(6,4,[1 2 3]);
ylim([0 1]);
plot(t,prey/max(prey))
hold on
plot(t,pred/max(pred))
%         plot(t,extmed/max(extmed+20));
xlim([0 t(length(t))])
text(-0.07 ,1.01,'a','Units','normalized','FontSize',13)
hold off 
xlabel('Time [d]','Interpreter','latex'); ylabel('Abundance','Interpreter','latex');

% plot prey wavelet
s21= subplot(6,4,8);
plot(log(mean(powerX,2))/log(mean(max(powerX),2)), log2(period));
text(-0.2,1.1,'d','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]); xlim([0 1.1]);       
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
s22 = subplot(6,4,[5 6 7]);
pcolor(t,log2(period),log2(powerX));
colormap parula
axis tight
shading flat
hold on
plot(t, log2(coi),'k')
ylim([min(log2(period)) max(log2(period))]);
text(-0.07,1.1,'c','Units','normalized','FontSize',13)
hold off
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse') 
ylabel('Period [d]','Interpreter','latex');

% plot predator wavelet spectrum
s31=subplot(6,4,12);
plot(log(mean(powerY,2))/log(mean(max(powerY),2)), log2(period));
text(-0.2,1.1,'f','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]); xlim([0 1.1]);
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
s32 = subplot(6,4,[9 10 11]);
pcolor(t,log2(period),log10(powerY));
colormap parula
axis tight
shading flat
hold on
plot(t, log2(coi),'k')
text(-0.07,1.1,'e','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]);
hold off
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
ylabel('Period [d]','Interpreter','latex');

% plot cross wavlet spectrum
s41= subplot(6,4,16);
plot(log(mean(powerXY,2))/log(mean(max(powerXY),2)), log2(period));
text(-0.2,1.1,'h','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]); xlim([0 1.1]);
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
xlabel('Power','Interpreter','latex')
s42= subplot(6,4,[13 14 15]);
pcolor(t,log2(period),log10(abs(powerXY)));
colormap parula
axis tight
shading flat
hold on
plot(t, log2(coi),'k')
text(-0.07,1.1,'g','Units','normalized','FontSize',13)
ylim([min(log2(period)) max(log2(period))]);
hold off
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
xlabel('Time [d]','Interpreter','latex')
ylabel('Period [d]','Interpreter','latex');

%plot coherence of wavelets
high = period(high)*ones(1,length(t));
low = period(low)*ones(1,length(t));
s51=subplot(6,4,20);
plot(mean(waveCo,2), log2(period));
text(-0.2,1.1,'j','Units','normalized','FontSize',13)
xlim([0 1]);
ylim([min(log2(period)) max(log2(period))]);
xlabel('Coherence','Interpreter','latex');
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
s52 = subplot(6,4,[17 18 19]);
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
text(-0.07 ,1.1,'i','Units','normalized','FontSize',13)
hold off
ylim([min(log2(period)) max(log2(period))]);
set(gca, 'YTickLabel', [2^2 2^3 4^2]);
set(gca, 'YDir','reverse')
xlabel('Time [d]','Interpreter','latex'); ylabel('Period [d]','Interpreter','latex');

%plot phase relationships
n = size(t,2);
s62= subplot(6,4,[21 22 23]);
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
plot(t,anglePlotXY/pi)
hold on
plot(t,anglePlotX3/pi)
plot(t,anglePlotX4/pi)
plot(t,anglePlotX5/pi)
text(-0.07,1.1,'k','Units','normalized','FontSize',13)
hold off
%         legend({'pred','egg-ratio','eggs','dead'},'Location','southeast')
xlim([0 t(length(t))])
xlabel('Time [d]','Interpreter','latex'); ylabel('$\phi$ $[\pi]$','Interpreter','latex')
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
s61del =  subplot(6,4,24);
s61 = polaraxes('Units',s61del.Units,'Position',s61del.Position);
delete(s61del);
hold on
polarhistogram(anglePlotXYi,9,'Normalization','pdf')
polarhistogram(anglePlotX3i,9,'Normalization','pdf')
polarhistogram(anglePlotX4i,9,'Normalization','pdf')
polarhistogram(anglePlotX5i,9,'Normalization','pdf')
text(-0.35,1.3,'l','Units','normalized','FontSize',13)
hold off

s11.Position = s11.Position + [0.03 0.04 -0.055 -0.02];
s12.Position = s12.Position + [0 0.04 0 -0.02];
s21.Position = s21.Position + [0.01 0.035 -0.03 -0.01];
s22.Position = s22.Position + [0 0.035 0 -0.01];
s31.Position = s31.Position + [0.01 0.045 -0.03 -0.01];
s32.Position = s32.Position + [0 0.045 0 -0.01];
s41.Position = s41.Position + [0.01 0.055 -0.03 -0.01];
s42.Position = s42.Position + [0 0.055 0 -0.01];
s51.Position = s51.Position + [0.01 0.05 -0.03 -0.01];
s52.Position = s52.Position + [0 0.05 0 -0.01];
s61.Position = s61.Position + [0 0.02 -0.01 -0.01];
s62.Position = s62.Position + [0 0.04 0 -0.02];
s = ['Figures/ExpC',c,mode];
%         print(gcf,'-depsc','-painters',s);
export_fig(s, '-pdf', '-transparent');    