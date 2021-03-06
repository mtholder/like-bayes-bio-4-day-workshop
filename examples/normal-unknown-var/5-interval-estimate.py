#!/usr/bin/env python
'''This program 2 or 3 arguments:
     a filename of a file with the data (one datum per line), and
     a standard deviation
Reports summary statistics and an MLE for mu
If a third argument is provided, it is interpreted as a value for mu
    and the likelihood is calculated.
'''
from __future__ import print_function
from scipy import optimize
import math

std_dev = None

def calc_ln_likelihood(data, m):
    global std_dev # do NOT use global variables in general!
    datum_ln_likelihood_list = []
    # Pr(data | m) = Product Pr(datum_i | m)
    # Pr(datum_i | m ) = (2 pi std_dev^2)^(-1/2) exp(-(datum_i - m)^2/(2 std_dev^2))
    coefficient = 1.0/(std_dev * math.sqrt(2*math.pi))
    ln_c = math.log(coefficient)
    var_factor = -1.0/(2.0*std_dev*std_dev)
    ln_L = 0.0
    for datum in data:
        diff = datum - m
        datum_ln_likelihood = ln_c + var_factor*diff*diff
        datum_ln_likelihood_list.append(datum_ln_likelihood)
        ln_L += datum_ln_likelihood
    return ln_L


def single_var_maximize_ln_L(data, ln_L_func, low_param, high_param):
    '''Takes `data` (the data set)
    `initial` is a starting parameter value
    `ln_L_func` is a function that should take the (data, parameter_value)
    and return a log-likelihood.
    `low_param` and `high_param` should be 2 numbers (with high_param > low_param)
        that can be used in the call to find a set of points that bracket the
        optimum.
    '''
    def scipy_ln_likelihood(x):
        '''SciPy minimizes functions. We want to maximize the likelihood. This
        function adapts our ln_likelihood function to the minimization context
        by returning the negative log-likelihood.
    
        We use this function with SciPy's minimization routine (minimizing the
        negative log-likelihood will maximize the log-likelihood).
        '''
        ln_L = ln_L_func(data, x)
        if VERBOSE:
            sys.stderr.write('In wrapper around ln_L_func with param = {m} lnL = {l}\n'.format(m=x, l=ln_L))
        return -ln_L
    # Now that we have wrapped the scipy_ln_likelihood, we can find an 
    #       approximation to the solution
    mle = optimize.brent(scipy_ln_likelihood,
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
        return abs(ln_L - ln_L_cutoff)

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
        return abs(ln_L - ln_L_cutoff)

    ub = find_upper_bound(data, ln_L_func, mle, ln_L_cutoff)
    return optimize.fminbound(scipy_ln_likelihood_upper_conf,
                              mle,
                              ub,
                              xtol=1e-8,
                              full_output=False)

if __name__ == '__main__':
    import sys
    import os
    VERBOSE = 'VERBOSE' in os.environ
    #   Read the command line arguments and set
    #   data as the list of data points
    #   std_dev as the sigma parameter for the normal distribution
    try:
        filename = sys.argv[1]
        std_dev = float(sys.argv[2])
        assert std_dev > 0
        user_mu = float(sys.argv[3]) if len(sys.argv) == 4 else None
        if not os.path.exists(filename):
            raise RuntimeError('"{f}" does not exist'.format(f=filename))
        with open(filename, 'r') as inp:
            data = [float(i.strip()) for i in inp]
    except Exception as x:
        sys.stderr.write('An error interpreting the input occurred. ')
        sys.stderr.write(__doc__)
        raise

    ln_l_function = calc_ln_likelihood

    n = len(data)
    mean = sum(data)/n
    print('sample mean = {m}'.format(m=mean))
    std_err = std_dev/math.sqrt(n)
    z_crit = 1.96
    margin_of_error = z_crit*std_err
    print('Z-statistic based 95% CI: ({l}, {u})'.format(l=(mean - margin_of_error),
                                                        u=(mean + margin_of_error)))
    user_ln_like = None
    if user_mu is not None:
        mu = user_mu
        user_ln_like = ln_l_function(data, mu)
        print('ln[Pr(data | mu={m})] = {l}'.format(m=mu, l=user_ln_like))
    mle = single_var_maximize_ln_L(data,
                                   ln_L_func=ln_l_function,
                                   low_param=min(data),
                                   high_param=max(data))
    mu = mle
    ln_like = ln_l_function(data, mu)
    print('MLE = {m}'.format(m=mu))
    print('ln[Pr(data | mu={m})] = {l}'.format(m=mu, l=ln_like))
    if user_ln_like is not None:
        ln_l_diff = 2*(ln_like - user_ln_like)
        print('2*(difference in lnL) = {}'.format(ln_l_diff))
        if ln_l_diff > 3.84:
            print('user-supplied parameter value is significantly worse using the chi-squared based cutoff')
        else:
            print('user-supplied parameter value is NOT significantly worse using the chi-squared based cutoff')
    else:
        chi_sq_critical = 3.84
        cutoff_diff = chi_sq_critical/2.0
        cutoff = ln_like - cutoff_diff
        lower_mu = find_lower_confidence_bound(data,
                                               ln_L_func=ln_l_function,
                                               mle=mle,
                                               ln_L_cutoff=cutoff)
        lower_ln_l = ln_l_function(data, lower_mu)
        print('lower-CI ln[Pr(data | lower_mu={m})] = {l}'.format(m=lower_mu, l=lower_ln_l))
        upper_mu = find_upper_confidence_bound(data,
                                               ln_L_func=ln_l_function,
                                               mle=mle,
                                               ln_L_cutoff=cutoff)
        upper_ln_l = ln_l_function(data, upper_mu)
        print('upper-CI ln[Pr(data | upper_mu={m})] = {l}'.format(m=upper_mu, l=upper_ln_l))
        print('Approx. 95% CI: ({l}, {u})'.format(l=lower_mu, u=upper_mu))


