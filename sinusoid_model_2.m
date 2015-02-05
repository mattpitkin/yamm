function y = sinusoid_model_2(t, parnames, parvals)

% function y = sinusoid_model_2(t, parnames, parvals)
%
% This function creates a sinusiod given by the equation:
%     y = amp * sin(phi0 + 2*pi*sum_n=0^4((1/(n+1)!)*f_n*(t-t0)^(n+1)))
% The input parameters are:
%
% t - the times (in sec) at which y will be calculated
% parnames - a cell array containing the parameter names.  
%     These can be in any order, but must include the following
%         {'amp', 'phi0', 't0', 'phi1', 'phi2', 'phi3', 'phi4', 'phi5'}
%     where amp is the amplitude, phi0 is the initial phase (in radians),
%     t0 is the phase and frequency epoch, and phi(n) are the initial
%     phase derivatives in radians at epochs 
% parvals - a cell array containing the values of the parameters
%     given in parnames.  These must be in the same order as in
%     parnames
%
%--------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
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

nphi = 0;

phimax = 6;
phis = zeros(phimax,1);

% extract parameter values
for ii=1:nparams
    switch parnames{ii}
        case 'amp'
            amp = parvals{ii};
        case 'phi0'
            phis(1) = parvals{ii};
            if nphi < 1
                nphi = 1;
            end
        case 't0'
            t0 = parvals{ii};
        case 'phi1'
            phis(2) = parvals{ii};
            if nphi < 2
                nphi = 2;
            end
        case 'phi2'
            phis(3) = parvals{ii};
            if nphi < 3
                nphi = 3;
            end
        case 'phi3'
            phis(4) = parvals{ii};
            if nphi < 4
                nphi = 4;
            end
        case 'phi4'
            phis(5) = parvals{ii};
            if nphi < 5
                nphi = 5;
            end
        case 'phi5'
            phis(6) = parvals{ii};
            if nphi < 6
                nphi = 6;
            end
    end
end

% convert phis into frequency and frequency derivatives
phis(2:end) = (phis(2:end)-phis(1))/(2*pi);

tt0 = t - t0;

phase = zeros(length(t),1);

if nphi > 1 
    % set matrix for conversion between phis and fs
    dtmatrix = zeros(nphi-1);
    for i=1:nphi-1
        dtmatrix(i,:) = i*(tt0(end)/(nphi-1));
        for j=2:nphi-1
            dtmatrix(i,j) = dtmatrix(i,j)^j;
        end
    end
    
    % convert between phases and frequencies
    fs = phis(2:end)/(dtmatrix);
    
    for i=1:nphi-1
        tt = tt0.^i;
        phase = phase + mod(fs(i)*tt, 1);
    end
end

% calculate sinusiod
y = amp * sin(2*pi*mod(phase, 1.) + phis(1));

end
