function scaled = rescale_parameters(prior, params)

% scaled = rescale_parameters(prior, params)
%
% This function will do the reverse of scale_parameters.

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
        scaled(i) = params(i)*(p4 - p3) + p3;
    elseif strcmpi(priortype, 'gaussian')
        p3 = prior{i,3};
        scaled(i) = params(i)*p4 + p3;
    elseif strcmpi(priortype, 'lognormal')
        scaled(i) = params(i);
    elseif strcmpi(priortype, 'jeffreys')
        scaled(i) = 10^(params(i)*(log10(p4) - log10(p3)) + log10(p3));
    end
end
