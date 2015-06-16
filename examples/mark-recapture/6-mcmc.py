#!/usr/bin/env python 
'''
Uses MCMC to estimate the probability of different values of theta (the number
    of double-headed coins) given data from independent trials that entail
    flipping all n coins.

A uniform prior over all values of theta is assumed.

Example invocation 
    python coin_contamination.py 1000000 4 1 4 3
means run:
    - 1 million iterations of MCMC
    - with n=4 coins in each trial
    - starting at the state of theta=1 double-headed coin out of the 4
    - two observations: 4 heads and 3 heads
    
'''
from __future__ import print_function
import random
import math
import sys
import os

####################################################
# The following 3 functions vary from application to application...
####################################################
PROPOSAL_WINDOW_SIZE = 0.01
def propose_new_state(theta):
    '''Return a randomly proposed state and the 
    log(Hastings ratio).
    '''
    # propose a state from among the adjacent states
    u_neg_half_to_half = random.random() - .50
    diff = PROPOSAL_WINDOW_SIZE*u_neg_half_to_half
    proposed_state = theta + diff
    if proposed_state < 0.0:
        proposed_state = -proposed_state
    if proposed_state > 1.0:
        proposed_state = 1.0 - (1.0 - proposed_state)
        proposed = theta + 1
    return proposed_state, 0.0

def calc_ln_likelihood(data, death_prob):
    #sys.stderr.write('death_prob = {d} ...\n'.format(d=death_prob))
    if death_prob < 0.0 or death_prob > 1.0:
        return float('-inf')
    global tag_loss_prob_list, initial_n
    # Pr(n_i | n_{i-1} ) = [(1 - tag_loss_prob)(1 - death_prob)]^(n_i)[1 - (1 - tag_loss_prob)(1 - death_prob)]^(n_{i-1} - n_i)
    prev_n = initial_n
    tlp_list = list(tag_loss_prob_list) # make a copy so that we don't change the global!
    ln_L = 0.0
    #sys.exit(str((initial_n, tag_loss_prob_list, data)))
    for curr_n in data:
        n_exited = prev_n - curr_n
        assert n_exited >= 0
        if tlp_list:
            tag_loss_prob = tlp_list.pop(0)
        recovery_prob = (1.0 - tag_loss_prob)*(1.0 - death_prob)
        exit_prob = 1.0 - recovery_prob
        if curr_n > 0:
            ln_L += curr_n * math.log(recovery_prob)
        if n_exited > 0:
            ln_L += n_exited * math.log(exit_prob)
        prev_n = curr_n
    #sys.stderr.write('death_prob = {d} lnL = {l}\n'.format(d=death_prob, l=ln_L))
    return ln_L

def calc_ln_prior(theta):
    '''Uniform prior.'''
    return 0.0
####################################################
# Here is our command-line interface
####################################################
if __name__ == '__main__':
    import sys
    import os
    VERBOSE = 'VERBOSE' in os.environ
    #   Read the command line arguments and set
    #   data as the list of data points
    #   std_dev as the sigma parameter for the normal distribution
    try:
        filename = sys.argv[1]
        num_it = int(sys.argv[2])
        initial_n = int(sys.argv[3])
        assert num_it > 0
        tag_loss_prob_list = [float(i) for i in sys.argv[4:]]
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


####################################################
# The code below is basically the same for any MCMC
####################################################
# This is MCMC using the Metropolis algorithm, except
# we often have to write out each visited state rather
# than keep a count of them in `mcmc_samples`
def metropolis_hastings(data, start_state, mcmc_samples, sample_freq):
    state = start_state
    ln_likelihood = calc_ln_likelihood(data, state)
    ln_prior = calc_ln_prior(state)
    sys.stderr.write("Iter\tlnL\tlnPrior\ttheta\n")
    for i in xrange(num_it):
        proposed_state, ln_hastings_ratio = propose_new_state(state)
        proposed_ln_prior = calc_ln_prior(proposed_state)
        proposed_ln_likelihood = calc_ln_likelihood(data, proposed_state)
        
        ln_prior_ratio = proposed_ln_prior - ln_prior
        ln_likelihood_ratio = proposed_ln_likelihood - ln_likelihood

        ln_posterior_ratio = ln_prior_ratio + ln_likelihood_ratio

        ln_acceptance_ratio = ln_posterior_ratio + ln_hastings_ratio
        
        if math.log(random.random()) < ln_acceptance_ratio:
            # Accept the move, update the state and current state info
            state = proposed_state
            ln_likelihood = proposed_ln_likelihood
            ln_prior = proposed_ln_prior
        if (i % sample_freq) == 0:
            mcmc_samples.append(state)
            sys.stderr.write("{i}\t{l}\t{p}\t{s}\n".format(i=i,
                                                           l=ln_likelihood,
                                                           p=ln_prior,
                                                           s=state))

start_state = random.random()
mcmc_samples = []
metropolis_hastings(data, start_state, mcmc_samples, 100)
mcmc_samples.sort()
mean = sum(mcmc_samples)/len(mcmc_samples)
print('posterior mean = {m}'.format(m=mean))
lower_ci_index = int(0.025*len(mcmc_samples))
upper_ci_index = int(0.975*len(mcmc_samples))
lower_ci = mcmc_samples[lower_ci_index]
upper_ci = mcmc_samples[upper_ci_index]

print('95% credible interval [{l}, {u}]'.format(l=lower_ci, u=upper_ci))
