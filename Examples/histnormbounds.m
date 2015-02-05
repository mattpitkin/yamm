function [n, x] = histnormbounds(y, m, low, high)

% [n, x] = histnormbounds(y, m, low, high)
%
% Create a normalised historam (normalised using the trapezium rule to
% calculated the area). The low and high values are the boundaries of the 
% histogrammed data, and ensure that the histogram gets properly normalised
% if close to the edges of the boundaries. m is the number of histogram 
% bins or a vector of bins values.

% set default bounds if none are set (low = -inf and high = inf)
if ~exist('low', 'var')
    low = -inf;
end

if ~exist('high', 'var')
    high = inf;
end

% check that data is within bounds
if min(y) < low && max(y) > high
    error('Data is outside of bounds');
end

% get histrogram
if isscalar(m)
    if m < 2
        error('Need more points in histogram');
    end
    
    [n, x] = hist(y, m);
else
    if length(m) < 2
        error('Need more points in histogram');
    end
    
    n = hist(y, m);
    x = m;
end

botbinwidth = x(2)-x(1);
topbinwidth = x(end)-x(end-1); % just in case x isn't uniform

% if histogram points are not close to boundaries (i.e. within a bin of the
% boundaries) then add zeros to histrogram edges
if x(1)-botbinwidth > low
    n = [0, n];
    x = [x(1)-botbinwidth, x];
end

if x(end)+topbinwidth < high
    n = [n, 0];
    x = [x, x(end)+topbinwidth];
end

% if the histogram is closer to the boundary edge than the bin width then
% set a new bin on the boundary with a value linearly extrapolated from the
% gradiant of the adjacent points
if x(1)-botbinwidth < low
    dx = x(1) - low;
    x = [low, x];
    
    dn = n(2)-n(1);
    
    nbound = n(1) - (dn/botbinwidth)*dx;
    
    n = [nbound, n];
end

if x(end)+topbinwidth > high
    dx = high - x(end);
    x = [x, high];
    
    dn = n(end)-n(end-1);
    
    nbound = n(end) + (dn/topbinwidth)*dx;
    
    n = [n, nbound];
end
    

% normalise n by area
n = n/trapz(x, n);
