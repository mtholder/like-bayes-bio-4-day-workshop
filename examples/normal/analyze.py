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
import sys
import os
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


def calc_likelihood(m):
    global std_dev, data
    datum_likelihood_list = []
    # Pr(data | m) = Product Pr(datum_i | m)
    # Pr(datum_i | m ) = (2 pi std_dev^2)^(-1/2) exp(-(datum_i - m)^2/(2 std_dev^2))
    coefficient = 1.0/(std_dev * math.sqrt(2*math.pi))
    var_factor = -1.0/(2.0*std_dev*std_dev)
    likelihood = 1.0
    for datum in data:
        diff = datum - m
        datum_likelihood = coefficient * math.exp(var_factor*diff*diff)
        datum_likelihood_list.append(datum)
        likelihood *= datum_likelihood
    return likelihood

def naive_calc_ln_likelihood(m):
    return math.log(calc_likelihood(m))


if user_mu is not None:
    ln_like = naive_calc_ln_likelihood(user_mu)
    print('ln[Pr(data | mu={m})] = {l}'.format(m=user_mu, l=ln_like))
else:
    n = len(data)
    mean = sum(data)/n
    print('mean = {m}'.format(m=mean))


