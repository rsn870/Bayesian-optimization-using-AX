%% Simulated Annealing Algorithm to find Optimal Site-ing of Rain Guages

% Find w that minimizes Vr subject to ||w||=K0

clear variables, close all

global precipp sigmae T N rainmean Sr varr Er alpha w Nensemb Neff
global precipp_k rainmean_k Er_k w_k Sr_k

% Neff = effective number of rain guages active currently (all gauges for
% which wi = 1)

%% Input data
% Years to include in time-series
yearstart = 1901; yearend = 2000;

% Months to average
monthstart = 6; monthend = 9;

% Cardinality of w, i.e., total number of instruments available

K = 100;
objfn = 2; %objfn = 1 (bias); objfn = 2 (variance); objfn = 3 (MSE)

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
% w = ones(N,1); % whether a measurement is present
% beta = ones(N,1); % weight to observations
% beta(1:2000) = 2;
%
% beta = beta/sum(w.*beta);
alpha = 0.5; %probability of missing observations
niter = 10; %number of iterations of Simulated Annealing

w_k = zeros(N,1);
windx_k = randperm(N,K); %index of w where sensors are present.
eps_k=0.1; %fraction of perturbation for new proposal at every itertion
w_k(windx_k)=1; %set indices with 1 as 
w = w_k;
Neff = nnz(w);
% [varmodel_k,biasmodel_k] = getmodelstats(w_k,beta,alpha);

varstore = nan(niter,1);
biasstore = nan(niter,1);
msestore = nan(niter,1);
betavarstore = nan(niter,Neff);
betabiasstore = nan(niter,Neff);
betamsestore = nan(niter,Neff);

windx_store = nan(niter,K);
%% Optimization parameters
A = []; b = [];
lb = zeros(Neff,1); ub = [];
Aeq = w(windx_k)'; beq = 1;
beta0 = ones(Neff,1)/Neff; % uniform weights
%beta0(~windx_k)=0;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

% Parameters for simulated annealing
kT = 1;
redfac = 0.9;
burnIN = 4;
kT_red_iter = 10;

%% Optimize based on objfn
if objfn == 1
    %% Initial values to the w-optimization
    precipp_k = precipp(:,windx_k);
    latlist_k = latlist(windx_k);
    lonlist_k = lonlist(windx_k);
    rainmean_k = NaN(T,1);
    for i = 1:T
        rainmean_k(i,1) = sum(precipp_k(i,:)'.*cosd(latlist_k))/sum(cosd(latlist_k));
    end
    Er_k = mean(precipp_k,1)';
    betabias = fmincon(@getmodelbias_k,beta0,A,b,Aeq,beq,lb,ub,[],options);
    biasmin = getmodelbias_k(betabias);
    
    for ii=1:niter
        if mod(ii,10)==0 || ii<=burnIN
            fprintf('ii=%d \n',ii)
        end
        windx_k = proposal(windx_k,N,eps_k); % Sample a new proposal. 
        w_k = zeros(N,1);
        w_k(windx_k)=1;
        windx_store(ii,:)=windx_k;
        w = w_k;
        Aeq = w(windx_k)'; beq = 1;
        precipp_k = precipp(:,windx_k);
        latlist_k = latlist(windx_k);
        lonlist_k = lonlist(windx_k);
        rainmean_k = NaN(T,1);
        for i = 1:T
            rainmean_k(i,1) = sum(precipp_k(i,:)'.*cosd(latlist_k))/sum(cosd(latlist_k));
        end
        Er_k = mean(precipp_k,1)';
        % Min. bias
        betabias_k = fmincon(@getmodelbias_k,betabias,A,b,Aeq,beq,lb,ub,[],options);
        biasmin_k = getmodelbias_k(betabias_k);
        
        dE = biasmin_k-biasmin;
        pdE = min(1,exp(-dE/kT));
        rE = rand;
        if rE<pdE || ii<burnIN
            biasstore(ii) = biasmin_k;
            betabiasstore(ii,:) = betabias_k;
            betabias = betabias_k;
            biasmin = biasmin_k;
        else
            biasstore(ii) = biasmin;
            betabiasstore(ii,:) = betabias;
        end
        if ii==burnIN
            kT = mean(biasstore(1:burnIN));
        elseif mod(ii,kT_red_iter)==0
            kT = redfac*kT;
        end
            
    end
end

if objfn == 2
    precipp_k = precipp(:,windx_k);
    latlist_k = latlist(windx_k);
    lonlist_k = lonlist(windx_k);
    rainmean_k = NaN(T,1);
    for i = 1:T
        rainmean_k(i,1) = sum(precipp_k(i,:)'.*cosd(latlist_k))/sum(cosd(latlist_k));
    end
    Er_k = mean(precipp_k,1)';
    Sr_k = Sr(windx_k,windx_k);
    betavar = fmincon(@getmodelvar_k,beta0,A,b,Aeq,beq,lb,ub,[],options);
    varmin = getmodelvar_k(betavar);
    
    for ii=1:niter
        if mod(ii,1)==0
            fprintf('ii=%d \n',ii)
        end
        windx_k = proposal(windx_k,N,eps_k);
        w_k = zeros(N,1);
        w_k(windx_k)=1;
        windx_store(ii,:)=windx_k;
        w = w_k;
        Aeq = w(windx_k)'; beq = 1;
        precipp_k = precipp(:,windx_k);
        latlist_k = latlist(windx_k);
        lonlist_k = lonlist(windx_k);
        rainmean_k = NaN(T,1);
        for i = 1:T
            rainmean_k(i,1) = sum(precipp_k(i,:)'.*cosd(latlist_k))/sum(cosd(latlist_k));
        end
        Er_k = mean(precipp_k,1)';
        Sr_k = Sr(windx_k,windx_k);
        
        % Min. variance
        betavar_k = fmincon(@getmodelvar_k,beta0,A,b,Aeq,beq,lb,ub,[],options);
        varmin_k = getmodelvar_k(betavar_k);
        
        dE = varmin_k-varmin;
        pdE = min(1,exp(-dE/kT));
        rE = rand;
        if rE<pdE || ii<burnIN
            varstore(ii) = varmin_k;
            betavarstore(ii,:) = betavar_k;
            betavar = betavar_k;
            varmin = varmin_k;
        else
            varstore(ii) = varmin;
            betavarstore(ii,:) = betavar;
        end
        
        if ii==burnIN
            kT = mean(varstore(1:burnIN));
        elseif mod(ii,kT_red_iter)==0
            kT = redfac*kT;
        end
        
        
    end
end

if objfn == 3
    betamse = fmincon(@getmodelmse_k,beta0,A,b,Aeq,beq,lb,ub,[],options);
    msemin = getmodelmse_k(betamse);
    for ii=1:niter
        if mod(ii,1)==0
            fprintf('ii=%d \n',ii)
        end
        windx_k = proposal(windx_k,N,eps_k);
        w_k = zeros(N,1);
        w_k(windx_k)=1;
        windx_store(ii,:)=windx_k;
        w = w_k;
        Aeq = w(windx_k)'; beq = 1;
        precipp_k = precipp(:,windx_k);
        latlist_k = latlist(windx_k);
        lonlist_k = lonlist(windx_k);
        rainmean_k = NaN(T,1);
        for i = 1:T
            rainmean_k(i,1) = sum(precipp_k(i,:)'.*cosd(latlist_k))/sum(cosd(latlist_k));
        end
        Er_k = mean(precipp_k,1)';
        
        % Min. MSE
        betamse_k = fmincon(@getmodelmse_k,beta0,A,b,Aeq,beq,lb,ub,[],options);
        msemin_k = getmodelmse_k(betamse_k);
        
        dE = msemin_k-msemin;
        pdE = min(1,exp(-dE/kT));
        rE = rand;
        if rE<pdE || ii<burnIN
            msestore(ii) = msemin_k;
            betamsestore(ii,:) = betamse_k;
            betamse = betamse_k;
            msemin = msemin_k;
        else
            msestore(ii) = msemin;
            betamsestore(ii,:) = betamse;
        end
        
        if ii==burnIN
            kT = mean(msestore(1:burnIN));
        elseif mod(ii,kT_red_iter)==0
            kT = redfac*kT;
        end
        
    end
end
