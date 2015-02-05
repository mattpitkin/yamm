function scaled = scale_parameters(prior, params)

% scaled = scale_parameters(prior, params)
%
% This function will scale parameters based on their priors. If a prior is
% uniform over a range then the parameter will be scaled, such that the
% range covers 0->1. If a prior is Gaussian then the parameter will be
% scaled such that it will be a Gaussian with zero mean and unit variance.
% If the prior is a function then the values will not be re-scaled. Also,
% at the moment (until I figure out how it should be done) lognormal prior
% parameters are not scaled.

lp = length(params);

scaled = zeros(lp,1);

for i=1:lp
    priortype = prior{i,2};
    p4 = prior{i,4};
    
    if strcmpi(priortype, 'function')
        scaled = params.';
        break;
    elseif strcmpi(priortype, 'uniform')
        p3 = prior{i,3};
        scaled(i) = (params(i) - p3)/(p4 - p3);
    elseif strcmpi(priortype, 'gaussian') || strcmpi(priortype, 'lognormal')
        p3 = prior{i,3};
        scaled(i) = (params(i) - p3)/p4;
    elseif strcmpi(priortype, 'lognormal')
        scaled(i) = params(i);
    elseif strcmpi(priortype, 'jeffreys')
        scaled(i) = (log10(params(i)) - log10(p3))/(log10(p4) - log10(p3));
    end
end
