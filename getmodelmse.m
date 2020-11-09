function [msemodel] = getmodelmse(beta)

global precipp sigmae T N rainmean Sr Er alpha w

%% Weights
gamma = w.*beta; % these are the effective weights used in the calculation

%% Variance

varsigma2 = (1/T*sum(precipp.*precipp,1) - alpha*(1/T*sum(precipp,1)).^2)';
Srnd = Sr;  I = eye(N); idiag = find(I); Srnd(idiag) = 0;

varm1 = 1/alpha*sum(gamma.^2.*varsigma2);
varm2 = gamma'*Srnd*gamma;
varm3 = 1/alpha*sigmae^2*(gamma'*gamma);
varm4 = -2*(1-alpha)/alpha*(gamma'*Er)*sum(gamma.^2.*Er);
varm5 = (1-alpha)/alpha*(gamma'*Er)^2*(gamma'*gamma);

varmodel = varm1 + varm2 + varm3 + varm4 + varm5; % equation (24) 

%% Bias
bias2 = (1-alpha)/alpha*((gamma'*Er)*(gamma'*gamma) - sum(gamma.^2.*Er));

% biassum = 0;
% for timecount = 1:T,
%     precippi = squeeze(precipp(timecount,:))';
%     bias1i = gamma'*precippi - rainmean(timecount);
%     biassum = biassum + (bias1i + bias2)^2; % equation (28)
% end
% 
% biasmodel = biassum / T;

a1 = gamma.*precipp';
a2 = sum(a1,1)';
bias1 = a2-rainmean;
biasmodel = mean((bias1+bias2).^2);

%% MSE
msemodel = biasmodel + varmodel; 
