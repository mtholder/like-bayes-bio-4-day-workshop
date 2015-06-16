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
# The following 3 functions calculate the likelihood for the strange
#   coin problem
####################################################
def n_choose_k(n, k):
    num, denom = 1, 1
    if 2*k > n:
        return n_choose_k(n, n-k)
    for i in range(k):
        num *= (n - i)
        denom *= (i + 1)
    nck = num/denom
    return nck

def fair_coin_ln_prob(h, n):
    return math.log(n_choose_k(n, h)) + n*math.log(0.5)

# since the range of each datum is small and the set 
#   of parameter values is small, we'll calculate 
#   the likelihood for every combination.
_cached_ln_likelihood_terms = []
def _calculate_all_ln_likelihood_factors(max_obs, num_coins):
    for n in range(max_obs + 1):
        p = []
        for theta in range(num_coins + 1):
            if theta > n:
                p.append(float('-inf'))
            else:
                p.append(fair_coin_ln_prob(n - theta, num_coins - theta))
        _cached_ln_likelihood_terms.append(p)


####################################################
# The following 3 functions vary from application to application...
####################################################
def propose_new_state(theta):
    '''Return a randomly proposed state and the 
    log(Hastings ratio).
    '''
    # propose a state from among the adjacent states
    if random.random() < 0.5:
        proposed = theta + 1
        if theta > num_coins:
            proposed = 0
    else:
        proposed = theta - 1
        if theta < 0:
            proposed = num_coins
    # our proposal is symmetric, so the Hastings ratio is 1
    # and the ln[Hastings ratio] = 0.0
    return proposed, 0.0

def calc_ln_likelihood(data, theta):
    '''ln_likelihood just as in the ML programs'''
    ln_like = 0.0
    for h in data:
        ln_like += _cached_ln_likelihood_terms[h][theta]
    return ln_like

def calc_ln_prior(theta):
    '''Uniform prior.'''
    return math.log(1.0/num_coins)
####################################################
# Here is our command-line interface
####################################################
try:
    num_it = int(sys.argv[1])
    assert(num_it > 0)
    num_coins = int(sys.argv[2])
    assert(num_coins > 0)
    state = int(sys.argv[3])
    assert(num_coins > 0)
    data = tuple([int(i) for i in sys.argv[4:]])
    assert(len(data) > 0)
    max_obs = max(data)
    assert(max_obs <= num_coins)
    assert(min(data) >= 0)
    # precalculate a lot of parts of the likelihood...
    _calculate_all_ln_likelihood_factors(max_obs, num_coins)
    # Whenever you write a program with a random number usage, it is wise to:
    #   1. report the pseudorandom number generator's seed, and
    #   2. supply some way of specifying the seed.
    # Here we use the "environment" that the script is run in. This is cryptic
    #   (not listed as a command line argument), but unobtrusive.  And we don't
    #   expect to use this option often in this script
    if 'MCMC_SEED' in os.environ:
        seed = float(os.environ['MCMC_SEED'])
    else:
        import time
        seed = time.time()
    random.seed(seed)
    msg_template = 'pseudorandom number seed = {s} (use the env. var. MCMC_SEED to set the seed).\n'
    sys.stderr.write(msg_template.format(s=repr(seed)))
except:
    sys.exit('Expecting:\npython coin_contamination.py <# iter> <# coins> <start state> <datum #1> <datum #2>...\n')


####################################################
# The code below is basically the same for any MCMC
####################################################
# This is MCMC using the Metropolis algorithm, except
# we often have to write out each visited state rather
# than keep a count of them in `mcmc_samples`
def metropolis_hastings(data, start_state, mcmc_samples):
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

        # Hastings ratio is 1.0, so we can ignore it...
        hastings_ratio = (0.5)/(0.5)

        ln_acceptance_ratio = ln_posterior_ratio + ln_hastings_ratio
        
        if math.log(random.random()) < ln_acceptance_ratio:
            # Accept the move, update the state and current state info
            state = proposed_state
            ln_likelihood = proposed_ln_likelihood
            ln_prior = proposed_ln_prior
        mcmc_samples[state] += 1
        sys.stderr.write("{i}\t{l}\t{p}\t{s}\n".format(i=i,
                                                       l=ln_likelihood,
                                                       p=ln_prior,
                                                       s=state))



# If there is a small number of discrete states we can just keep track
#    the number of iterations that we spend in each state
mcmc_samples = [0]*(num_coins + 1)
metropolis_hastings(data, state, mcmc_samples)
    
            
print("Posterior probabilities from MCMC")
for state in range(num_coins + 1):
    print(state, float(mcmc_samples[state])/num_it)

print("\nTrue Posterior probabilities (calculated analytically)")
likelihood_list = [math.exp(calc_ln_likelihood(data, i)) for i in range(num_coins + 1)]
marginal_prob = sum(likelihood_list)
for state, likelihood in enumerate(likelihood_list):
    print(state, likelihood/marginal_prob)
print("\nTransition Probabilities:\n                       From")
print("     " + "    ".join(["%7d" % i for i in range(num_coins + 1)]))
for num_dh in range(num_coins + 1):
    ind_below = num_dh - 1
    ind_above = num_dh + 1 if num_dh < num_coins else 0
    same_state = 0.0
    likelihood = likelihood_list[num_dh]
    if likelihood == 0:
        ti_prob_above, ti_prob_below = 0.5, 0.5
    else:
        like_below, like_above = likelihood_list[ind_below], likelihood_list[ind_above]
        if like_above > likelihood:
            ti_prob_above = 0.5
        else:
            ti_prob_above = 0.5*like_above/likelihood
            same_state += (0.5 - ti_prob_above)
        if like_below > likelihood:
            ti_prob_below = 0.5
        else:
            ti_prob_below = 0.5*like_below/likelihood
            same_state += (0.5 - ti_prob_below)
    ti_probs = [0.0] * (num_coins + 1)
    ti_probs[num_dh] = same_state
    ti_probs[ind_below] = ti_prob_below
    ti_probs[ind_above] = ti_prob_above
    print("%-4d  %s" % (num_dh, " ".join([" %7f " % d for d in ti_probs])))
