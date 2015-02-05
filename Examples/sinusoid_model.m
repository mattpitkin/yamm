function y = sinusoid_model(t, parnames, parvals)

% function y = sinusoid_model(t, parnames, parvals)
%
% This function creates a sinusiod given by the equation:
%     y = amp * sin(phi0 + 2*pi*sum_n=0^4((1/(n+1)!)*f_n*(t-t0)^(n+1)))
% The input parameters are:
%
% t - the times (in sec) at which y will be calculated
% parnames - a cell array containing the parameter names.  
%     These can be in any order, but must include the following
%         {'amp', 'phi0', 't0', 'f0', 'f1', 'f2', 'f3', 'f4'}
%     where amp is the amplitude, phi0 is the initial phase (in radians),
%     t0 is the phase and frequency epoch, and f(n) are the initial
%     frequency and its time derivative (in Hz/s^n)
% parvals - a cell array containing the values of the parameters
%     given in parnames.  These must be in the same order as in
%     parnames
%
%--------------------------------------------------------------------------
%           This is the format required by mcmc_sampler.m.
%--------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

nparams = lpn;

% extract parameter values
for ii=1:nparams
    switch parnames{ii}
        case 'amp'
            amp = parvals{ii};
        case 'phi0'
            phi0 = parvals{ii};
        case 't0'
            t0 = parvals{ii};
        case 'f0'
            f0 = parvals{ii};
        case 'f1'
            f1 = parvals{ii};
        case 'f2'
            f2 = parvals{ii};
        case 'f3'
            f3 = parvals{ii};
        case 'f4'
            f4 = parvals{ii}; 
    end
end

fs = [f0, f1, f2, f3, f4];

tt0 = t - t0;
phase = zeros(length(t),1);
for i=1:5
    tt = tt0.^i;
    phase = phase + mod(fs(i)*tt/factorial(i), 1);
end

% calculate sinusiod
y = amp * sin(2*pi*mod(phase, 1.) + phi0);

end
