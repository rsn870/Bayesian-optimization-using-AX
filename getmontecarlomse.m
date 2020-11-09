function [msedat] = getmontecarlomse(beta)

global precipp sigmae T N rainmean w alpha Nensemb

%% Weights
gamma = w.*beta; % these are the effective weights used in the calculation

%% Monte Carlo Simulation
varensemblist = NaN(Nensemb,1);

for outcount = 1:Nensemb,
    
    pmeanoptm = NaN(T,1);
    pmeanoptmne = NaN(T,1);
    
    if rem(outcount,500) == 0,
        simnum = outcount
    end
    
    for count = 1:T,
        scount = binornd(1,alpha,N,1);
        precippi = (precipp(count,:) + randn(1,N)*sigmae)';
        pmeanoptm(count) = sum(gamma.*precippi.*scount)/sum(gamma.*scount);
        
        precippine = precipp(count,:)';
        pmeanoptmne(count) = sum(gamma.*precippine.*scount)/sum(gamma.*scount);
    end
    
    if outcount == 1,
        pmeanoptensemb = pmeanoptm;
    else
        pmeanoptensemb = pmeanoptensemb + pmeanoptm;
    end
  
    varensemblist(outcount,1) = var(pmeanoptm);
    
end

pmeanoptensemb = pmeanoptensemb / Nensemb;

vardat = mean(varensemblist);
biasdat = mean((pmeanoptensemb-rainmean).^2);

msedat = vardat + biasdat;