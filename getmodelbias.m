function [biasmodel] = getmodelbias(beta)

% global precipp T rainmean Er alpha w
global precipp rainmean Er alpha w 
%% Weights
gamma = w.*beta; % these are the effective weights used in the calculation

%% Bias
bias2 = (1-alpha)/alpha*((gamma'*Er)*(gamma'*gamma) - sum(gamma.^2.*Er));
%%
% t2 = tic;
% biassum = 0;
% for timecount = 1:T
%     precippi = squeeze(precipp(timecount,:))';
%     bias1i = gamma'*precippi - rainmean(timecount);
%     biassum = biassum + (bias1i + bias2)^2; % equation (28)
% end
% 
% biasmodel_old = biassum / T;
% toc(t2);
% %%
% keyboard
% %%
% t1 = tic;
% a1 = bsxfun(@(x,y) x.*y,gamma,precipp');
% a2 = sum(a1,1)';
% bias1 = a2-rainmean;
% biasmodel = mean((bias1+bias2).^2);
% toc(t1)
%%

%%
% t3 = tic;
a1 = gamma.*precipp';
a2 = sum(a1,1);
bias1 = a2'-rainmean;
biasmodel = mean((bias1+bias2).^2);
% toc(t3)
%%
% keyboard
% %%
% t4 = tic;
% a4a = gpuArray(precipp);
% a1 = gamma.*a4a';
% a2 = sum(a1,1);
% bias1 = a2'-rainmean;
% biasmodel = mean((bias1+bias2).^2);
% toc(t4)