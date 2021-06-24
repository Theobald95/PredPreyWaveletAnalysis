%   Code straight forwardly translated to MATLAB from 1d Wavelet transform from Torrence and Compo

%   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
%
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made. This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of the above software in any publications:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Research Systems, Inc.                 Climate Diagnostics Center
%  4990 Pearl East Circle                 325 Broadway R/CDC1
%  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
%  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
%----------------------------------------------------------------------------
%%
function [wave, period, scale, coi] = wavelet(y,dt,pad,dj,s0,j1,mother,param)

    n1 = length(y);
    if s0 == -1; s0=2*dt; end
    if dj == -1; dj = 1. / 4.; end
    if j1 == -1; j1=round((log(n1*dt/s0)/log(2))/dj); end
    if mother == -1; mother = "MORLET"; end

    %....construct time series to analyze, pad if necessary
    x = y - mean(y);
    if pad == 1
        x = TaperData(x);
        x = PaddData(x);
    end
    n = length(x);

    %....construct wavenumber array used in transform [Eqn(5)]
    k = 1:round(n/2);  
    k = k .* ((2 * pi)/(n*dt));
    k = [0., k, -k(round((n-1)/2):-1:1)];
    %
    %....compute FFT of the (padded) time series
    f = fft(x);    % [Eqn(3)]

    %....construct SCALE array & empty PERIOD & WAVE arrays
    scale = s0*2 .^ ((0:j1)*dj);
    wave = zeros(j1+1,n);   % define the wavelet array

    % loop through all scales and compute transform
    for a1 = 1:j1+1
       [daughter,fourier_factor,coi,~] = wave_bases(mother,k,scale(a1),param);
       wave(a1,:) = ifft(f.*daughter(1:length(daughter)-1));  % wavelet transform[Eqn(4)]
    end

    period = fourier_factor*scale;
    %coi = coi*dt*[1e-5; collect(1:((n1+1)/2-1)); collect((n1/2-1):-1:1); 1e-5]  % COI [Sec.3g]
    coi = coi*dt*[1e-5, (1:((n1+1)/2-1)), ((n1/2-1):-1:1), 1e-5];  % COI [Sec.3g]
    wave = wave(:,1:n1);  % get rid of padding before returning
    
end
%%
function [daughter, fourier_factor, coi, dofmin] = wave_bases(mother,k,scale,param)
mother = upper(mother);
n = length(k);
if mother == "MORLET"  %-----------------------------------  Morlet
   if param == -1; param = 6.; end
   k0 = param;
   expnt = -(scale*k - k0).^2 / 2 .* heaviside(k);
   norm = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);     % total energy=N   [Eqn(7)]
   daughter = norm*exp(expnt);
   daughter = daughter .* heaviside(k);                  % Heaviside step function
   fourier_factor = (4.0*pi)/(k0 + sqrt(2.0 + k0^2)); % Scale-->Fourier [Sec.3h]
   coi = fourier_factor/sqrt(2.0);                    % Cone-of-influence [Sec.3g]
   dofmin = 2;                                        % Degrees of freedom
 elseif mother == "PAUL"  %-----------------------------------  Paul
   if param == -1; param = 4. ; end
   m = param;
   expnt = -(scale.*k).* heavi.(k);
   norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
   daughter = norm*((scale.*k).^m).*exp(expnt);
   daughter = daughter.* heaviside(k);        % Heaviside step function
   fourier_factor = 4*pi/(2.0*m+1);
   coi = fourier_factor*sqrt(2.0);
   dofmin = 2;
 elseif mother == "DOG"  %-----------------------------------  DOG
   if param == -1; param = 2.; end
   m = param;
   expnt = -(scale.*k).^2 ./ 2.0;
   norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
   daughter = -norm*(im^m)*((scale.*k).^m).*exp(expnt);
   fourier_factor = 2*pi*sqrt(2 ./ (2*m+1));
   coi = fourier_factor/sqrt(2);
   dofmin = 1;
 else
   error("Mother must be one of MORLET,PAUL,DOG")
end
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
barlett = 1:N;
barlett = 1-abs((barlett-0.5*N)/(0.5*N));
tapData = zeros(size(X));
for i = 1:size(X,1)
    tapData(i,:)= X(i,:).*barlett;
end
tapData = (N/ sum(barlett))*tapData;
end