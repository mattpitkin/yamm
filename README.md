# yamm
Yet Another Matlab MCMC code

Here is another Matlab (and potentially Octave compatible) code for performing Markov chain Monte Carlo parameter estimation.

The code can use a variety of proposal functions including the "stretch" and "walk" affine invariant ensemble samplers of [Goodman & Weare](http://msp.org/camcos/2010/5-1/p04.xhtml). It also allows a variety of prior functions (or a user defined function) for the parameters.

The distribution includes some examples of setting up likelihood functions and models, which I will document at some point!

For Matlab code doing a very similar job see e.g. [The MCMC Hammer](http://www.mathworks.com/matlabcentral/fileexchange/49537-the-mcmc-hammer---affine-invariant-mcmc-sampler) by Aslak Grinsted, or look elsewhere on the [Matlab File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/) or [Google](https://www.google.co.uk/#q=matlab+mcmc)! Or, alternatively if you're more into Python (as I am too!) check out [emcee](http://dan.iel.fm/emcee/current/) (as described in Forman-Mackay et al., [arXiv:1202.3665](http://arxiv.org/abs/1202.3665)).

For similar Matlab-based Bayesian parameter estimation and evidence evaluation you may also want to checkout the [Matlab implementation of MultiNest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest/frs/).

&copy; Matthew Pitkin, 2015
