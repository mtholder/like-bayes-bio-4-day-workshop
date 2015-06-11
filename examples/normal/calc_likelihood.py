#!/usr/bin/env python
'''This program 2 or 3 arguments:
     a filename of a file with the data (one datum per line), and
     a standard deviation
Reports summary statistics and an MLE for mu
If a third argument is provided, it is interpreted as a value for mu
    and the likelihood is calculated.
'''
from __future__ import print_function
import math

std_dev = None

def calc_likelihood(data, m):
    global std_dev # do NOT use global variables in general!
    datum_likelihood_list = []
    # Pr(data | m) = Product Pr(datum_i | m)
    # Pr(datum_i | m ) = (2 pi std_dev^2)^(-1/2) exp(-(datum_i - m)^2/(2 std_dev^2))
    coefficient = 1.0/(std_dev * math.sqrt(2*math.pi))
    var_factor = -1.0/(2.0*std_dev*std_dev)
    likelihood = 1.0
    for datum in data:
        diff = datum - m
        datum_likelihood = coefficient * math.exp(var_factor*diff*diff)
        datum_likelihood_list.append(datum_likelihood)
        likelihood *= datum_likelihood
    return likelihood

def naive_calc_ln_likelihood(data, m):
    return math.log(calc_likelihood(data, m))

if __name__ == '__main__':
    import sys
    import os
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

    ln_l_function = naive_calc_ln_likelihood

    n = len(data)
    mean = sum(data)/n
    print('mean = {m}'.format(m=mean))
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
    else:
        sys.exit('This version of the program requires a third argument')
