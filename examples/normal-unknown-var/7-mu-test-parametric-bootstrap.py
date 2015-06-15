#!/usr/bin/env python
'''This program  1 or 2 arguments:
     a filename of a file with the data (one datum per line), and
It reports summary statistics and an MLE for mu
If a second argument is provided, it is interpreted as a value for mu
    and the likelihood is calculated.
'''
from __future__ import print_function
from scipy import optimize
import random
import math

def calc_ln_likelihood(data, mu, std_dev):
    if VERBOSE:
        sys.stderr.write('calc_ln_likelihood(data, mu={m}, sd={s}) = ...\n'.format(m=mu, s=std_dev))
    # Pr(data | mu) = Product Pr(datum_i | mu)
    # Pr(datum_i | mu ) = (2 pi std_dev^2)^(-1/2) exp(-(datum_i - mu)^2/(2 std_dev^2))
    coefficient = 1.0/(std_dev * math.sqrt(2*math.pi))
    ln_c = math.log(coefficient)
    var_factor = -1.0/(2.0*std_dev*std_dev)
    ln_L = 0.0
    for datum in data:
        diff = datum - mu
        datum_ln_likelihood = ln_c + var_factor*diff*diff
        ln_L += datum_ln_likelihood
    if VERBOSE:
        sys.stderr.write('calc_ln_likelihood(data, mu={m}, sd={s}) = {l}\n'.format(m=mu, s=std_dev, l=ln_L))
    return ln_L


def maximize_ln_L_for_fixed_mu(data, ln_L_func, mu):
    '''Takes `data` (the data set)
    `ln_L_func` is a function that should take the (data, parameter_value)
    and return a log-likelihood.
    `mu` is the second arg to `ln_L_func` it is fixed for this optimization.
    `low_param` and `high_param` should be 2 numbers (with high_param > low_param)
        that can be used in the call to find a set of points that bracket the
        optimum.
    '''
    def sigma_scipy_ln_likelihood(sd):
        '''SciPy minimizes functions. We want to maximize the likelihood. This
        function adapts our ln_likelihood function to the minimization context
        by returning the negative log-likelihood.
    
        We use this function with SciPy's minimization routine (minimizing the
        negative log-likelihood will maximize the log-likelihood).
        '''
        if sd <= 0.0:
            return float('inf')
        ln_L = ln_L_func(data, mu, sd)
        if VERBOSE:
            sys.stderr.write('In wrapper around ln_L_func with mu = {m} lnL = {l}\n'.format(m=mu, l=ln_L))
        return -ln_L
    diff_list = []
    for datum in data:
        diff_list.append(abs(mu - datum))
    diff_list.sort()
    smallest_diff = diff_list[0]
    largest_diff = diff_list[-1]
    # Now that we have wrapped the scipy_ln_likelihood, we can find an 
    #       approximation to the solution
    mle = optimize.brent(sigma_scipy_ln_likelihood,
                         brack=(smallest_diff, largest_diff),
                         tol=1e-8,
                         full_output=False)
    return mle

def optimize_sigma_ln_L(data, x):
    s = maximize_ln_L_for_fixed_mu(data, ln_L_func=calc_ln_likelihood, mu=x)
    return calc_ln_likelihood(data, x, s)

def maximize_mu_ln_L_over_all_sigma(data, ln_L_func, low_param, high_param):
    '''Takes `data` (the data set)
    `ln_L_func` is a function that should take the (data, parameter_value)
    and return a log-likelihood.
    `low_param` and `high_param` should be 2 numbers (with high_param > low_param)
        that can be used in the call to find a set of points that bracket the
        optimum.
    '''
    def mu_scipy_ln_likelihood(x):
        '''Here we perform an optimization to find the MLE of sigma
        for our current value of mu'''
        # a variance is sort of an average of the squared differences
        #   between the mean that the variates.
        # So a good guess at a bracketing range is the smallest absolvute value
        #   of a difference and the largest absolute value of the difference
        sigma_mle = maximize_ln_L_for_fixed_mu(data, ln_L_func, x)
        ln_L = ln_L_func(data, x, sigma_mle)
        if VERBOSE:
            sys.stderr.write('for mu={m} nuisance sigma_mle is {s} lnL = {l}\n'.format(m=x, s=sigma_mle, l=ln_L))
        return -ln_L
    # Now that we have wrapped the scipy_ln_likelihood, we can find an 
    #       approximation to the solution
    mle = optimize.brent(mu_scipy_ln_likelihood,
                         brack=(low_param, high_param),
                         tol=1e-8,
                         full_output=False)
    return mle

def find_lower_bound(data,
                     ln_L_func,
                     mle,
                     ln_L_cutoff):
    eps = 0.1 # to avoid crashes if the MLE= 0
    lb = mle *.95 - eps
    lb_ln_L = ln_L_func(data, lb)
    if VERBOSE:
        sys.stderr.write('In lower conf finder lb = {m} lnL = {l}\n'.format(m=lb, l=lb_ln_L))
    while lb_ln_L > ln_L_cutoff:
        value_diff = mle - lb
        lb = lb - value_diff
        lb_ln_L = ln_L_func(data, lb)
        if VERBOSE:
            sys.stderr.write('In lower conf finder lb = {m} lnL = {l}\n'.format(m=lb, l=lb_ln_L))
    return lb

def find_upper_bound(data,
                     ln_L_func,
                     mle,
                     ln_L_cutoff):
    eps = 0.1 # to avoid crashes if the MLE= 0
    ub = mle *1.05 + eps
    ub_ln_L = ln_L_func(data, ub)
    if VERBOSE:
        sys.stderr.write('In upper conf finder ub = {m} lnL = {l}\n'.format(m=ub, l=ub_ln_L))
    while ub_ln_L > ln_L_cutoff:
        value_diff = ub - mle
        ub = ub + value_diff
        ub_ln_L = ln_L_func(data, ub)
        if VERBOSE:
            sys.stderr.write('In upper conf finder ub = {m} lnL = {l}\n'.format(m=ub, l=ub_ln_L))
    return ub

def find_lower_confidence_bound(data,
                                ln_L_func,
                                mle,
                                ln_L_cutoff):
    def scipy_ln_likelihood_lower_conf(x):
        if x > mle:
            return float('inf')
        ln_L = ln_L_func(data, x)
        if VERBOSE:
            sys.stderr.write('In lower conf wrapper around ln_L_func with param = {m} lnL = {l}\n'.format(m=x, l=ln_L))
        return (ln_L - ln_L_cutoff)**2

    lb = find_lower_bound(data, ln_L_func, mle, ln_L_cutoff)
    return optimize.fminbound(scipy_ln_likelihood_lower_conf,
                              lb,
                              mle,
                              xtol=1e-8,
                              full_output=False)

def find_upper_confidence_bound(data,
                                ln_L_func,
                                mle,
                                ln_L_cutoff):
    def scipy_ln_likelihood_upper_conf(x):
        if x < mle:
            return float('inf')
        ln_L = ln_L_func(data, x)
        if VERBOSE:
            sys.stderr.write('In upper conf wrapper around ln_L_func with param = {m} lnL = {l}\n'.format(m=x, l=ln_L))
        return (ln_L - ln_L_cutoff)**2

    ub = find_upper_bound(data, ln_L_func, mle, ln_L_cutoff)
    return optimize.fminbound(scipy_ln_likelihood_upper_conf,
                              mle,
                              ub,
                              xtol=1e-8,
                              full_output=False)


def calculate_lrt(data, null_mu):
    mle = maximize_mu_ln_L_over_all_sigma(data,
                                          ln_L_func=calc_ln_likelihood,
                                          low_param=min(data),
                                          high_param=max(data))
    mu = mle
    sigma_mle = maximize_ln_L_for_fixed_mu(data,
                                          ln_L_func=calc_ln_likelihood,
                                          mu=mu)
    ln_like = calc_ln_likelihood(data, mu, sigma_mle)
    null_sigma_mle = maximize_ln_L_for_fixed_mu(data,
                                                ln_L_func=calc_ln_likelihood,
                                                mu=null_mu)
    null_ln_like = calc_ln_likelihood(data, null_mu, null_sigma_mle)
    mle = {'mu': mu,
           'sigma': sigma_mle,
           'lnL': ln_like,
          }
    null = {'mu': null_mu,
           'sigma': null_sigma_mle,
           'lnL': null_ln_like,
          }
    lrt = ln_like - null_ln_like
    return mle, null, lrt


def simulate(n, mu, sigma):
    sim_data = []
    for i in xrange(n):
        sim_data.append(random.normalvariate(mu, sigma))
    return sim_data

def parametric_bootstrap(n_reps, n, null, lrt):
    null_mu = null['mu']
    num_more_extreme = 0
    for i in xrange(n_reps):
        sys.stderr.write('    sim rep {r}\n'.format(r=i))
        sim_data = simulate(n, null['mu'], null['sigma'])
        sim_mle, sim_null, sim_lrt = calculate_lrt(sim_data, null_mu)
        if sim_lrt > lrt:
            num_more_extreme += 1
    return num_more_extreme/float(n_reps)

if __name__ == '__main__':
    import sys
    import os
    VERBOSE = 'VERBOSE' in os.environ
    try:
        filename = sys.argv[1]
        user_mu = float(sys.argv[2])
        if not os.path.exists(filename):
            raise RuntimeError('"{f}" does not exist'.format(f=filename))
        with open(filename, 'r') as inp:
            data = [float(i.strip()) for i in inp]
    except Exception as x:
        sys.stderr.write('An error interpreting the input occurred. ')
        sys.stderr.write(__doc__)
        raise
    user_ln_like = None
    mle, null, lrt = calculate_lrt(data,
                                   null_mu=user_mu)
    mu_mle = mle['mu']
    sigma_mle = mle['sigma']
    max_ln_L = mle['lnL']
    print('MLE = mu={m}, sigma={s}'.format(m=mu_mle, s=sigma_mle))
    print('ln[Pr(data | mu={m}, sigma={s})] = {l}'.format(m=mu_mle, s=sigma_mle, l=max_ln_L))
    mu_null = null['mu']
    sigma_null = null['sigma']
    null_ln_L = null['lnL']
    print('ln[Pr(data | mu_null={m}, sigma={s})] = {l}'.format(m=mu_null, s=sigma_null, l=null_ln_L))
    print('LRT = {l}'.format(l=lrt))
    approx_p = parametric_bootstrap(1000, len(data), null, lrt)
    print('P is approximately = {p}'.format(p=approx_p))
    



