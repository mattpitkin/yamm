function [sampleOut, logPOut, logPriorOut, ar] = draw_mcmc_sample(sampleIn, cholmat, ...
    logPIn, logPriorIn, prior, data, likelihood, model, parnames, extraparvals, ...
    temperature, dimupdate)

% function [sampleOut, logPOut, logPriorOut, ar] = draw_mcmc_sample(sampleIn, cholmat, ...
%    logPIn, logPriorIn, prior, data, likelihood, model, Nmcmc, parnames, ...
%    extraparvals, temperature, dimupdate)
%
% If sampleIn contains one sample then this function will draw a new
% multi-dimensional sample using the proposal distribution defined by
% cholmat. If sampleIn contains multiple samples (an ensmeble) then it will
% use the stretch method (see Goodman & Weare, CAMCOS, Vol. 5 (2010), 1,
% 65–80 or Foreman-Mackay et al, http://arxiv.org/abs/1202.3665) to update
% the ensemble of "walkers".
% 
% In the first case the MCMC will use a Students-t (with N=2 degrees of
% freedom) proposal distribution based on the Cholesky decomposed
% covariance matrix, cholmat. Whilst in the second case the stretch
% parameter, z, will be drawn from a distribution between 1/3 and 3.
%
% The posterior probability will be calculated using the prior and the
% likelihood. extraparvals is a vector of additional parameters needed by
% the model. The sample to be returned will be based on the
% Metropolis-Hastings rejection criterion. In the single sample case the
% proposal used will always be symmetric, so the proposal ratio is unity.
% Whilst in the ensmeble case the proposal ratio with be given by z^N,
% where N is the number of parameters/dimensions.
% 
% The dimupdate specifies the number of dimensions/parameters that you wish
% to update with each newly drawn sample. If this is 0 or greater than the
% actual number of dimensions then all parameters will be updated,
% otherwise the number of specified parameters will be updated, with the
% updated parameters chosen at random.
%
% The returned array ar contains one(s) if the sample(s) are accepted,
% or 0 if rejected.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l2p = 0.5*log(2*pi); % useful constant
Ndegs = 2; % degrees of freedom of Students't distribution

ss = size(sampleIn);

Nens = ss(1);  % number of samples in ensemble
Npars = ss(2); % number of parameters in a sample

ar = zeros(Nens,1);

sampletmp = zeros(Nens, Npars);
sampleOut = zeros(Nens, Npars);
logPOut = zeros(Nens,1);
logPriorOut = zeros(Nens,1);
propratio = zeros(Nens,1);

% if not using an ensemble of samples do this
if Nens == 1
    % draw new points from mulitvariate Gaussian distribution 
    gasdevs = randn(Npars, 1);
    sampletmp(1,:) = (cholmat*gasdevs)';

    % check how many of the dimensions you want to update
    if dimupdate > 0 && dimupdate < Npars
        r = randperm(Npars);
        sampletmpnew = zeros(1, Npars);
        sampletmpnew(1, r(1:dimupdate)) = sampletmp(1, r(1:dimupdate));
        sampletmp = sampletmpnew;
    end

    % calculate chi-square distributed value
    chi = sum(randn(Ndegs,1).^2);
            
    % add value onto old sample
    sampletmp = sampleIn + sampletmp*sqrt(Ndegs/chi);
else
    % update the ensemble of samples with stretch or walk moves (use
    % stretch 75% of the time and walk the other 25%)
    
    if rand < 0.75 % do stretch moves
        % draw stretch value from between 1/a and a (uniformly in log(a))
        maxscale = 3;
        R = rand(Nens,1);
        lm = log(maxscale);
        logscale = 2.*lm*R - lm;
        scale = exp(logscale); % scale factor
        propratio = logscale*Npars; % proposal ratio
        
        if dimupdate > 0 && dimupdate < Npars
            r = randperm(Npars);
            propratio = logscale*dimupdate; % re-calculate proposal ratio
        end
    
        % loop through ensemble and update parameters
        for j=1:Nens
            % pick another random sample from the ensemble to us in the update
            i = randi(Nens);
            while i == j
                i = randi(Nens);
            end
        
            sampletmp(j,:) = scale(j)*(sampleIn(j,:)-sampleIn(i,:));
        
            % check how many of the dimensions you want to update
            if dimupdate > 0 && dimupdate < Npars
                sampletmpnew = zeros(1, Npars);
                sampletmpnew(1, r(1:dimupdate)) = sampletmp(1, r(1:dimupdate));
                sampletmp(j,:) = sampletmpnew;
            end
        
            % add value on to samples
            sampletmp(j,:) = sampleIn(i,:) + sampletmp(j,:);
        end
    else % do walk moves
        sample_size = 3; % number of samples to pick from ensemble
        
        propratio = zeros(Nens, 1); % symmetric proposals
        
        % loop through ensemble and update parameters
        for j=1:Nens
            % choose random samples from ensemble
            idxs = randperm(Nens, sample_size);
            
            % get centre of mass of points
            com = mean(sampleIn(idxs,:));
            
            % generate univariate random Gaussian variables for each
            % parameteter
            uv = randn(sample_size,1);
            
            % get steps
            sampletmp(j,:) = (sampleIn(idxs,:)-repmat(com,sample_size,1))'*uv;
            
            
            % check how many of the dimensions you want to update
            if dimupdate > 0 && dimupdate < Npars
                r = randperm(Npars);
                sampletmpnew = zeros(1, Npars);
                sampletmpnew(1, r(1:dimupdate)) = sampletmp(1, r(1:dimupdate));
                sampletmp(j,:) = sampletmpnew;
            end
            
            % update sample
            sampletmp(j,:) = sampleIn(j,:) + sampletmp(j,:);
        end
    end
end

for i=1:Nens
    % check sample(s) is within the (scaled) prior
    newPrior = 0;
    outofbounds = 0;

    for j=1:Npars
        priortype = prior{j,2};
        p3 = prior{j,3};
        p4 = prior{j,4};
        
        funcprior = 0;
        
        % check if prior is a function
        if strcmpi(priortype, 'function')
            priorfunc = prior{j,5};
            if isa(priorfunc, 'function_handle')
                funcprior = 1;
            elseif exist(priortype, 'file') == 2
                funcprior = 1;
                priorfunc = str2func(priorfunc);
            else
                error('Prior function is not defined properly');
            end
        elseif strcmpi(priortype, 'uniform')
            behaviour = prior{j,5};
                
            dp = 1;
                
            if sampletmp(i,j) < 0 || sampletmp(i,j) > 1
                if strcmpi(behaviour, 'reflect')
                    % reflect the value from the boundary
                    sampletmp(i,j) = 1 - mod(sampletmp(i,j), dp);
                elseif strcmpi(behaviour, 'cyclic')
                    % wrap parameter from one side to the other
                    while sampletmp(i,j) > 1
                        sampletmp(i,j) = sampletmp(i,j) - 1;
                    end
                        
                    while sampletmp(i,j) < 0
                        sampletmp(i,j) = sampletmp(i,j) + 1;
                    end
                else
                    outofbounds = 1;
                    break;
                end
            end
                
            pv = -log(p4-p3);
            newPrior = newPrior + pv;
            
        elseif strcmpi(priortype, 'gaussian')
            pv = -l2p - sampletmp(i,j)^2/2;
            newPrior = newPrior + pv;
        elseif strcmpi(priortype, 'lognormal')
            pv = log(lognpdf(sampletmp(i,j), p3, p4));
            newPrior = newPrior + pv;
        elseif strcmpi(priortype, 'jeffreys')
            behaviour = prior{j,5};

            dp = 1;
               
            if sampletmp(i,j) < 0 || sampletmp(i,j) > 1
                if strcmp(behaviour, 'reflect')
                    % reflect the value from the boundary
                    sampletmp(i,j) = 1 - mod(sampletmp(i,j), dp);
                else
                    outofbounds = 1;
                    break;
                end
            end
                
            pv = -log(10^(sampletmp(i,j)*(log10(p4) - log10(p3)) + log10(p3)));
            newPrior = newPrior + pv;
        end
        
        if funcprior
            sc = rescale_parameters(prior, sampletmp(i,:));
            pv = feval(priorfunc, parnames, cat(1, num2cell(sc), extraparvals));
            newPrior = newPrior + pv;
            
            % reject point if prior is zero
            if newPrior == -inf
                outofbounds = 1;
            end
            
            break; % prior function should be prior for all parameters so break out of loop
        end
    end
    
    if outofbounds % reject point
        sampleOut(i,:) = sampleIn(i,:);
        logPOut(i) = logPIn(i);
        logPriorOut(i) = logPriorIn(i);
    else
        % rescale sample back to its proper range for likelihood
        sc = rescale_parameters(prior, sampletmp(i,:));
    
        % get likelihood of new point
        logLnew = feval(likelihood, data, model, parnames, cat(1, num2cell(sc), ...
            extraparvals));
    
        % get posterior probability
        logPnew = logLnew + newPrior;
        
        ratio = ((logPnew - logPIn(i)) * temperature + propratio(i));
        X = log(rand);
        
        % accpet/reject point using Metropolis-Hastings criterion
        if X > ratio
            % reject the new sample
            sampleOut(i,:) = sampleIn(i,:);
            logPOut(i) = logPIn(i);
            logPriorOut(i) = logPriorIn(i);
        else
            % accept the new sample
            sampleOut(i,:) = sampletmp(i,:);
            logPOut(i) = logPnew;
            logPriorOut(i) = newPrior;
            ar(i) = 1;
        end
    end
end

return
