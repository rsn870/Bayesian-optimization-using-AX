%% For 0.25 x 0.25 deg field
%% Calculate statistics (bias, variance, and MSE) 
%% Given weights beta and w
clear all, close all

global precipp sigmae T N rainmean Sr varr Er alpha w Nensemb

%% Input data
% Years to include in time-series
yearstart = 1901; yearend = 2000;

% Months to average
monthstart = 6; monthend = 9;

%% Input parameters
% measurement standard deviation
sigmae = 1; % mm / day

% size of Monte carlo ensemble
Nensemb = 5000;

%% load & read data, 0.25 deg x 0.25 deg
load raindat0_25_deg % lonlist latlist hrrainmat

% hrrainmat is {year, day, location}

%% Consider a subset of location
% This is for testing purposes only
% Comment out for paper
Nuse = 4964;
lonlist = lonlist(1:Nuse);
latlist = latlist(1:Nuse);
hrrainmat = hrrainmat(:,:,1:Nuse);


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
for i = 1:T,
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
%beta = ones(N,1); % weight to observations

%beta = beta/sum(w.*beta);

%% Loop for different values of alpha
%alphalist = linspace(0.2,1,17); % probability of missing observation
alphalist = 0.5; % probability of missing observation

varmat = NaN(numel(alphalist),2); biasmat = NaN(numel(alphalist),2);

for i = 1:numel(alphalist),
    alpha = alphalist(i)
    
    %% Optimization parameters
    A = []; b = [];
    lb = zeros(N,1); ub = [];
    Aeq = w'; beq = 1;
    beta0 = ones(N,1)/sum(w); % uniform weights
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    
    %% Optimization from model
    % Min. bias
%     betabias = fmincon(@getmodelbias,beta0,A,b,Aeq,beq,lb,ub,[],options); 
%     biasmin = getmodelbias(betabias);
    
    % Min. variance
    betavar = fmincon(@getmodelvar,beta0,A,b,Aeq,beq,lb,ub,[],options); 
    varmin = getmodelvar(betavar);
    
    % Min. MSE
    betamse = fmincon(@getmodelmse,beta0,A,b,Aeq,beq,lb,ub,[],options); 
    msemin = getmodelmse(betamse);
    
    minbiaslist(i) = biasmin;
    minvarlist(i) = varmin;
    minmselist(i) = msemin;
    
    %% Optimization from Monte Carlo simulation
    
    mcflag = 0;
    
    if mcflag == 1
        
        % Min. bias
        %betabiasmc = fmincon(@getmontecarlobias,beta0,A,b,Aeq,beq,lb,ub,[],options);
        biasminmc = getmontecarlobias(betabias);
        
        % Min. variance
        %betavarmc = fmincon(@getmontecarlovar,beta0,A,b,Aeq,beq,lb,ub,[],options);
        varminmc = getmontecarlovar(betavar);
        
        % Min. MSE
        %betamsemc = fmincon(@getmontecarlomse,beta0,A,b,Aeq,beq,lb,ub,[],options);
        mseminmc = getmontecarlomse(betamse);
        
        minbiaslistmc(i) = biasminmc;
        minvarlistmc(i) = varminmc;
        minmselistmc(i) = mseminmc;
        
    end
end