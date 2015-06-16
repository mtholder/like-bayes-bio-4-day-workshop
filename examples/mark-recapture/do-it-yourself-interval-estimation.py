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

EPSILON = 1e-6
tag_loss_prob_list = None

def calc_ln_likelihood(data, death_prob):
    global tag_loss_prob_list, initial_n
    # Pr(n_i | n_{i-1} ) = [(1 - tag_loss_prob)(1 - death_prob)]^(n_i)[1 - (1 - tag_loss_prob)(1 - death_prob)]^(n_{i-1} - n_i)
    prev_n = initial_n
    tlp_list = list(tag_loss_prob_list) # make a copy so that we don't change the global!
    ln_L = 0.0
    for curr_n in data:
        n_exited = prev_n - curr_n
        if tlp_list:
            tag_loss_prob = tlp_list.pop(0)
        ###########################################
        # DO NOT EDIT ABOVE HERE
        # curr_n is the number of tags seen in this month.
        # n_exited is the difference between this number and the previous month.
        # tag_loss_prob is the probability of a tag that is working at the 
        #   beginning of the previous month falls off or stops working during
        #   the previous month
        # death_prob is the probability that an individual alive at the 
        #   beginning of the previous month is dead at this point in time.
        ###########################################
        #print('death', death_prob)
        prob_recovery = (1.0 - death_prob)*(1.0 - tag_loss_prob)
        #print('prob_recovery', prob_recovery)
        prob_not_rec = 1.0 - prob_recovery
        ln_L_this_datum = curr_n*math.log(prob_recovery) + n_exited * math.log(prob_not_rec)
        ###########################################
        # DO NOT EDIT BELOW HERE
        ###########################################
        ln_L += ln_L_this_datum
        prev_n = curr_n
    print('death = ', death_prob)
    #assert False
    return ln_L

# ALL of the functions below here are VERY similar to the 
#   code in the "normal" example

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
    # Here we use a bounded optimization, becaue the death prob can't be 1 or
    #       less than 0.0
    mle = optimize.fminbound(scipy_ln_likelihood,
                             0.0,
                             1.0 - EPSILON,
                             xtol=1e-8,
                             full_output=False)
    return mle


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

    lb = 0.0
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

    ub = 1.0 - EPSILON
    return optimize.fminbound(scipy_ln_likelihood_upper_conf,
                              mle,
                              ub,
                              xtol=1e-8,
                              full_output=False)

# The argument handling has to change some from the 
#   "normal" example because the user supplies different
#   information. 
# The only new concept here is storing a series of tag
#   loss probabilities in a list so that we can let the 
#   user tell us a variable number of probabilities.
if __name__ == '__main__':
    import sys
    import os
    VERBOSE = 'VERBOSE' in os.environ
    #   Read the command line arguments and set
    #   data as the list of data points
    #   std_dev as the sigma parameter for the normal distribution
    try:
        filename = sys.argv[1]
        initial_n = float(sys.argv[2])
        assert initial_n > 0
        tag_loss_prob_list = [float(i) for i in sys.argv[3:]]
        for tlp in tag_loss_prob_list:
            assert tlp >= 0.0
            assert tlp <= 1.0
        if not os.path.exists(filename):
            raise RuntimeError('"{f}" does not exist'.format(f=filename))
        with open(filename, 'r') as inp:
            data = [float(i.strip()) for i in inp]
    except Exception as x:
        sys.stderr.write('An error interpreting the input occurred. ')
        sys.stderr.write(__doc__)
        raise

    ln_l_function = calc_ln_likelihood

    mle = single_var_maximize_ln_L(data,
                                   ln_L_func=ln_l_function,
                                   low_param=min(data),
                                   high_param=max(data))
    mu = mle
    ln_like = ln_l_function(data, mu)
    print('MLE = {m}'.format(m=mu))
    print('ln[Pr(data | mu={m})] = {l}'.format(m=mu, l=ln_like))
    if True:
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


