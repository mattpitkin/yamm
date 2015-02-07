# yamm
## Yet Another Matlab MCMC code

Here is another Matlab (and potentially Octave compatible) code for performing Markov chain Monte Carlo parameter estimation. Other MCMC codes are available.

The code can use a variety of proposal functions including the "stretch" and "walk" affine invariant ensemble samplers of [Goodman & Weare](http://msp.org/camcos/2010/5-1/p04.xhtml). It also allows a variety of prior functions (or a user defined function) for the parameters.

## Setup

You can use the code by simply adding the ``src`` directory to your Matlab path, e.g. in ``bash`` in GNU/Linux just type

    export MATLABPATH=${MATLABPATH}:path_to_repository/src

Alternatively run the provided ``setup.py`` script with

    . ./setup.sh
    
to add the paths to ``MATLABPATH``, or use the ``addpath`` command in Matlab, or just do it through the window options.

If trying this with Octave just set the ``OCTAVEPATH`` environment variable instead, or again use the ``addpath`` command.

## Example

Here I'll give a couple of examples (well one at the moment) of setting up the likelihood and model functions and the parameters for estimation. The style of the likelihood and model functions given here are how they should be set up in general. Other examples are given in the ``Examples`` directory.

### Sinusoid example

An example of running the code for a sinusoidal model and Gaussian likelihood is given here. First let's define the model function. Note that the model function should take three variables:
 1. a variable, in this case ``t`` (a vector of time stamps), that could be a single value, vector, or matrix, giving the points at which the model is evaluated;
 2. a cell array (``parnames``) containing the names of any variables needed to define the model, where the names are recognised by the ``switch`` statement in the function;
 3. a cell array (``parvals``) containing the values of each of the variables in ``parnames``, which much be included in the same order. These variables can be anything - floats, integers, strings, vectors, matrices - as long as the model function knows what to do with them. However, any variables that are being estimated must just be individual real numbers.

```matlab
function y = sinusoid_model(t, parnames, parvals)

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

% extract parameter values
for i=1:lpn
    switch parnames{i}
        case 'amp' % sinusoid amplitude
            amp = parvals{i};
        case 't0' % the phase epoch
            t0 = parvals{i}
        case 'phi0' % sinusoid initial phase at t0
            phi0 = parvals{i};
        case 'freq' % sinusoid frequency (Hz)
            freq = parvals{i};
    end
end

% generate the sinusoid at the values t
y = amp * sin(2*pi*freq*(t-t0) + phi0);

```

Now let's define the log likelihood function, which we'll take as a Gaussian with a stationary noise variance. Note that the log likelihood function should take in four variables:
 1. a ``data`` variable containing any 'data' information. In this case ``data`` is a cell array with the first cell containing a vector of time stamps of the data samples, the second cell containing the actual data samples, and the third cell containing the data variance, but data could contain any other information required by the likelihood function, provided the likelihood function knows what to do with it;
 2. a ``model`` function (either a function handle or function name string) e.g. the ``sinusoid_model`` above;
 3. a cell array ``parnames`` containing the model parameter names to be passed to the model (as above);
 4. a cell array ``parvals`` containing the model parameter values to be passed to the model (as above).

```matlab
function logL = logL_gaussian(data, model, parnames, parvals)

% check whether model is a string or function handle
if ischar(model)
    fmodel = str2func(model);
elseif isa(model, 'function_handle')
    fmodel = model;
else
    error('Error... Expecting a model function!');
end

% get the data from the data cell array
t = data{1}; % vector of time stamps at which the data is given
y = data{2}; % vector of data (e.g. noise + signal)
C = data{3}; % the value of the noise variance

N = length(t);
invC = 1/C;
lDetC = N*log(C);

% evaluate the model
md = feval(fmodel, x, parnames, parvals);

% if the model returns a NaN then set the likelihood to be zero (e.g. loglikelihood to be -inf)
if isnan(md)
    logL = -inf;
    return;
end

% calculate the log likelihood
logL = -0.5*(y - md)' * invC * (y - md);
logL = logL - 0.5*N*log(2*pi) - 0.5*lDetC;

if isnan(logL)
    error('Error: log likelihood is NaN!');
end

return

```

Now let's create a script that creates some simulated data (a sinusoid added to Gaussian noise) and then uses the MCMC to extract the posteriors of two of the sinusoid parameters (the amplitude and the phase) assuming the frequency and phase epoch are known.
```matlab
% SIMULATE DATA
% discrete times
dt = 1/5; % sampling rate of signal (Hz)
tlen = 6; % length of signal (seconds)
t = (0:dt:(tlen-dt))';

% inject a simulated sinusoidal signal model
amp = 10.0; % signal amplitude
phi0 = 2.3; % initial phase of signal (rads)
f0 = 0.788634; % signal frequency
t0 = 0; % phase epoch

% create signal injection
s = sinusoid_model(t, {'amp', 'phi0', 't0', 'freq'}, {amp, phi0, t0, freq});

% gaussian noise model (white)
sigma2 = 16; % variance of Gaussian noise
n = sqrt(sigma2) * randn(length(t), 1);

% data (t, y, covariance matrix)
y = s+n; % signal + noise
data{1} = t;
data{2} = y;
data{3} = sigma2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define nested sampling parameters
Nmcmc = 50000; % MCMC samples
Nburnin = 50000; % burn in

likelihood = @logL_gaussian;
model = @sinusoid_model;

% setup the prior (THIS IS IMPORTANT)
% we want to estimate the 'amp' and 'phi0' parameters, so these go in a prior variable
% that defines the type of prior and range of the prior - in this case uniform priors
% over a given minimum and maximum range
prior = {'amp', 'uniform', 0, 20, 'fixed'; ...
         'phi0', 'uniform', 0, 2*pi, 'cyclic'};
         
% we need to define the other parameters needed for the model in the extraparams 
% cell array variable
extraparams = {'freq', freq; ...
               't0', t0};

% RUN THE MCMC (using the ensemble sampler with 10 'walkers')
[post_samples, logP] = mcmc_sampler(data, likelihood, model, prior, ...
    extraparams, 'Nmcmc', Nmcmc, 'Nburnin', Nburnin, ...
    'NensembleStretch', 10);

```

## Alternatives

For Matlab code doing a very similar job see e.g. [The MCMC Hammer](http://www.mathworks.com/matlabcentral/fileexchange/49537-the-mcmc-hammer---affine-invariant-mcmc-sampler) by Aslak Grinsted, the [MCMC toolbox for Matlab](http://helios.fmi.fi/~lainema/mcmc/) by Marko Laine, or look elsewhere on the [Matlab File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/) or [Google](https://www.google.co.uk/#q=matlab+mcmc)! Or, alternatively if you're more into Python (as I am too!) check out [emcee](http://dan.iel.fm/emcee/current/) (as described in Forman-Mackay et al., [arXiv:1202.3665](http://arxiv.org/abs/1202.3665)).

For similar Matlab-based Bayesian parameter estimation and evidence evaluation you may also want to checkout the [Matlab implementation of MultiNest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest/frs/).

&copy; Matthew Pitkin, 2015
