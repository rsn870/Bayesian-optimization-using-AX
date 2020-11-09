function [varmodel] = getmodelvar_k(beta)

global precipp_k sigmae T Neff Sr_k Er_k alpha

%% Weights
gamma = beta; % these are the effective weights used in the calculation

%% Variance
% t1 = tic;
% varsigma2 = (1/T*sum(precipp.*precipp,1) - alpha*(1/T*sum(precipp,1)).^2)';
% Srnd = Sr;  I = logical(eye(N)); Srnd(I)=0;
% %idiag = find(I); Srnd(idiag) = 0;
% %gtmp = gamma'*gamma;
% 
% varm1 = 1/alpha*sum(gamma.^2.*varsigma2);
% varm2 = gamma'*Srnd*gamma;
% varm3 = 1/alpha*sigmae^2*(gamma'*gamma);
% varm4 = -2*(1-alpha)/alpha*(gamma'*Er)*sum(gamma.^2.*Er);
% varm5 = (1-alpha)/alpha*(gamma'*Er)^2*(gamma'*gamma);
% 
% varmodel = varm1 + varm2 + varm3 + varm4 + varm5; % equation (24) 
% toc(t1)
% %%
% t1 = tic;
% varsigma2 = (1/T*sum(precipp.*precipp,1) - alpha*(1/T*sum(precipp,1)).^2)';
% Srnd = Sr;  
% Srnd = Srnd - diag(diag(Srnd));
% %Srnd(1:1+N:end)=0;
% toc(t1)
% %I = logical(eye(N)); Srnd(I)=0;
% %idiag = find(I); Srnd(idiag) = 0;
% 
% gdot = gamma'*gamma;
% gEr = gamma'*Er;
% 
% t2 = tic;
% varm1 = 1/alpha*sum(gamma.^2.*varsigma2);
% varm2 = gamma'*Srnd*gamma;
% varm3 = 1/alpha*sigmae^2*(gdot);
% varm4 = -2*(1-alpha)/alpha*(gEr)*sum(gamma.^2.*Er);
% varm5 = (1-alpha)/alpha*(gEr)^2*(gdot);
% toc(t2)
% 
% varmodel = varm1 + varm2 + varm3 + varm4 + varm5; % equation (24) 
% toc(t1)
%%
% t1 = tic;
varsigma2 = (1/T*sum(precipp_k.*precipp_k,1) - alpha*(1/T*sum(precipp_k,1)).^2)';
Srnd = Sr_k;  
%Srnd = Srnd - diag(diag(Srnd));
%Srnd(1:1+N:end)=0;
% toc(t1)
%I = logical(eye(N)); Srnd(I)=0;
%idiag = find(I); Srnd(idiag) = 0;

gdot = gamma'*gamma;
gEr = gamma'*Er_k;

% t2 = tic;
varm1 = 1/alpha*sum(gamma.^2.*varsigma2);
varm2 = gamma'*Srnd*gamma-sum((gamma.^2).*diag(Srnd));
varm3 = 1/alpha*sigmae^2*(gdot);
varm4 = -2*(1-alpha)/alpha*(gEr)*sum(gamma.^2.*Er_k);
varm5 = (1-alpha)/alpha*(gEr)^2*(gdot);
% toc(t2)

varmodel = varm1 + varm2 + varm3 + varm4 + varm5; % equation (24) 
% toc(t1)

%%
% keyboard
