%% Simulated Annealing Algorithm to find Optimal Citing of Rain Guages

% Find w that minimizes Vr subject to ||w||=K0

clear all, close all

global precipp sigmae T N rainmean Sr varr Er alpha w Nensemb

%% Input data
% Years to include in time-series
yearstart = 1901; yearend = 2000;

% Months to average
monthstart = 6; monthend = 9;

% Cardinality of w, total number of instruments available

K = 200;

%% Input parameters
% measurement standard deviation
sigmae = 1; % mm / day

% size of Monte carlo ensemble
Nensemb = 5000;

%% load & read data, 0.25 deg x 0.25 deg
load raindat0_25_deg % lonlist latlist hrrainmat

% hrrainmat is {year, day, location}

%% calculate average
datmat = permute(hrrainmat,[2 1 3]); % {day, year, location}

%% compute indices
numdaysinmonth = [30 31 30 31 31 30 31 30]; % number of days in April through November
cumsumdays = cumsum(numdaysinmonth);
cumsumdays = cat(2,0,cumsumdays);

daystartindex = cumsumdays(monthstart-3) + 1;
dayendindex = cumsumdays(monthend-3 + 1);

%% average datamat over chosen months and include time-series for selected years
datmatuse = datmat(daystartindex:dayendindex,yearstart-1900:yearend-1900,:);

% average over first index, which is days in chosen months
datmatuse2d = squeeze(mean(datmatuse,1));

precipp = datmatuse2d; % this is the data matrix [t=year,location]

N = size(precipp,2); % number of locations
T = size(precipp,1); % number of years

%% cosine-weighted mean
rainmean = NaN(T,1);
for i = 1:T
    rainmean(i,1) = sum(precipp(i,:)'.*cosd(latlist))/sum(cosd(latlist));
end

% rainmean is cosine-latitude weighted mean rainfall rate over chosen months, in mm/day

%% Covariance matrix

% mean observation
pmean = mean(precipp,1); % this is the mean of observations across years

% covariance estimation
Xp = precipp - repmat(pmean,[T 1]); % anomaly matrix

Sv = 1/(T-1)*(Xp'*Xp); % field covariance matrix

Sr = Sv + sigmae^2*eye(N); % observation covariance matrix

varr = diag(Sr);

%% Temporal mean
Er = pmean';

%% Observation placement and weights
w = ones(N,1); % whether a measurement is present
% beta = ones(N,1); % weight to observations
% beta(1:2000) = 2;
%
% beta = beta/sum(w.*beta);
alpha = 0.5;
niter = 100;

w_k = zeros(N,1);
windx_k = randperm(N,K);
eps_k=0.1;
w_k(windx_k)=1;
w = w_k;
% [varmodel_k,biasmodel_k] = getmodelstats(w_k,beta,alpha);

varstore = nan(niter,1);
biasstore = nan(niter,1);
msestore = nan(niter,1);
betavarstore = nan(niter,N);
betabiasstore = nan(niter,N);
betamsestore = nan(niter,N);

windx_store = nan(niter,K);
%% Optimization parameters
A = []; b = [];
lb = zeros(N,1); ub = [];
Aeq = w'; beq = 1;
beta0 = ones(N,1)/sum(w); % uniform weights
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%% Initial values to the w-optimization



betamse = fmincon(@getmodelmse,beta0,A,b,Aeq,beq,lb,ub,[],options);
msemin = getmodelmse(betamse);

for ii=1:niter
    if mod(ii,1)==0
        fprintf('ii=%d \n',ii)
    end
    windx_k = proposal(windx_k,N,eps_k);
    w_k = zeros(N,1);
    w_k(windx_k)=1;
    windx_store(ii,:)=windx_k;
    w = w_k;
    Aeq = w';
    %     [varmodel_k,biasmodel_k] = getmodelstats(w_k,beta,alpha);
    % Min. bias
    betabias = fmincon(@getmodelbias,beta0,A,b,Aeq,beq,lb,ub,[],options);
    biasmin = getmodelbias(betabias);
    
    % Min. variance
    betavar = fmincon(@getmodelvar,beta0,A,b,Aeq,beq,lb,ub,[],options);
    varmin = getmodelvar(betavar);
    
    % Min. MSE
    betamse = fmincon(@getmodelmse,beta0,A,b,Aeq,beq,lb,ub,[],options);
    msemin = getmodelmse(betamse);
    
    
    varstore(ii) = varmin;
    biasstore(ii) = biasmin;
    msestore(ii) = msemin;
    
    betavarstore(ii,:) = betavar;
    betabiasstore(ii,:) = betabias;
    betamsestore(ii,:) = betamse;
    
    
end
