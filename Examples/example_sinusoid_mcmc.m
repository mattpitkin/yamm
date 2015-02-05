% example: estimate amplitude and initial phase of a 
% sinusoidal signal in additive white gaussian noise

global verbose
verbose = 1;

% SIMULATE DATA
% discrete times
dt = 1/5; % sampling rate of signal (Hz)
tlen = 6; % length of signal (seconds)
t = (0:dt:(tlen-dt))';

% inject a simulated sinusoidal signal model
amp = 10.0; % signal amplitude
phi0 = 2.3; % initial phase of signal (rads)
f0 = 0.788634; % signal frequency
f1 = 0;
f2 = 0;
f3 = 0;
f4 = 0;
t0 = 0;

% create signal injection
s = sinusoid_model(t, ...
    {'amp', 'phi0', 't0', 'f0', 'f1', 'f2', 'f3', 'f4'}, ...
    {amp, phi0, t0, f0, f1, f2, f3, f4});

% gaussian noise model (white)
sigma2 = 16; % variance of Gaussian noise
n = sqrt(sigma2) * randn(length(t), 1);

% data (t, y, covariance matrix)
y = s+n; 
data{1} = t;
data{2} = y;
data{3} = ones(1, length(t))*sigma2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define nested sampling parameters
Nmcmc = 50000; % MCMC samples
Nburnin = 50000; % burn in

likelihood = @logL_gaussian;
model = @sinusoid_model;
prior = {'amp', 'uniform', 0, 20, 'fixed'; 'phi0', 'uniform', 0, 2*pi, 'cyclic'};
extraparams = {'f0', f0; 't0', t0; 'f1', f1; 'f2', f2; 'f3', f3; 'f4', f4}; % fixed signal freqs

temperature = 1;

% run MCMC sampler
[post_samples, logP] = mcmc_sampler(data, likelihood, model, prior, ...
    extraparams, 'Nmcmc', Nmcmc, 'Nburnin', Nburnin, 'temperature', ...
    temperature, 'NensembleStretch', 10);

% plot the samples
[namp, xamp] = histnormbounds(post_samples(:,1), 20, prior{1,3}, prior{1,4});
figure;
plot(xamp, namp, 'b', [amp amp], [0, max(namp)+0.1*max(namp)], 'k--');
set(gca, 'fontsize', 14, 'fontname', 'helvectica');
xlabel('Amplitude', 'fontsize', 14, 'fontname', 'avantgarde');
ylabel('Probability Density', 'fontsize', 14, 'fontname', 'avantgarde');

fprintf(1, 'Mean and standard deviation amplitude = %.6f +/- %.6f\n', ...
    mean(post_samples(:,1)), std(post_samples(:,1)));

[nphi0, xphi0] = histnormbounds(post_samples(:,2), 20, prior{2,3}, prior{2,4});
figure;
plot(xphi0, nphi0, 'b', [phi0 phi0], [0, max(nphi0)+0.1*max(nphi0)], 'k--');
set(gca, 'fontsize', 14, 'fontname', 'helvectica');
xlabel('\phi_0 (rads)', 'fontsize', 14, 'fontname', 'avantgarde');
ylabel('Probability Density', 'fontsize', 14, 'fontname', 'avantgarde');

fprintf(1, 'Mean and standard deviation of phi0 = %.6f +/- %.6f (rads)\n', ...
    mean(post_samples(:,2)), std(post_samples(:,2)));
