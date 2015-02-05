function [Neff, acl] = effective_sample_size(series)

% get the auto-correlation of a given series
x = series - mean(series);

z = conv(x, x', 'full');

z = ifftshift(z);

acf = z(1:length(series));

acf = acf/acf(1);

% estimate the auto-correlation length
M = 5;
K = 2;

acf(2:end) = acf(2:end)*2;

imax = floor(length(acf)/K);

cacf = cumsum(acf);

s = (1:(length(cacf)+1))/M;

idxs = find(cacf(1:imax) < s(1:imax));

% check whether idx is empty
if isempty(idxs)
    acl = 1;
else
    acl = s(idxs(1));
end
 

% effective number of uncorrelated samples
Neff = floor(length(series)/acl);