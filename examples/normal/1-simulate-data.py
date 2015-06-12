#!/usr/bin/env python
'''Take 3 arguments:
    sample_size mean standard_deviation
and writes to standard output a simulated realization from the normal distribution
(one variate per line).
'''
# If you were trying to do this in one line, you could use:
#    print('\n'.join([str(random.normalvariate(mean, std_dev)) for i in xrange(n)]))

# In a command-line program, we have 2 output "streams" that we can write to
# You should used sys.stderr for status messages and error messages; and
#   use sys.stdout for the "real" output of the program. "print" writes to
#   sys.stdout
# Some notes on the output below:
#   Writing to the ouput stream requires that we convert our value to a "string"
#   We  can create a template string with content like {v} to indicate where a
#     a value for v should go in the string form.
#   The \n represents a newline.
#   The format method of a string treats the string as a template and performs
#     variable substitution on the string using the values provided in the format
#     function call.
from __future__ import print_function
import random
import sys
import os
try:
    arguments = sys.argv[1:]
    n = int(arguments[0])
    assert n > 0
    mean = float(arguments[1])
    std_dev = float(arguments[2])
    assert std_dev >= 0
    # Whenever you write a program with a random number usage, it is wise to:
    #   1. report the pseudorandom number generator's seed, and
    #   2. supply some way of specifying the seed.
    # Here we use the "environment" that the script is run in. This is cryptic
    #   (not listed as a command line argument), but unobtrusive.  And we don't
    #   expect to use this option often in this script
    if 'SIMULATE_DATA_SET_SEED' in os.environ:
        seed = float(os.environ['SIMULATE_DATA_SET_SEED'])
    else:
        import time
        seed = time.time()
    random.seed(seed)
    msg_template = 'pseudorandom number seed = {s} (use the env. var. SIMULATE_DATA_SET_SEED to set the seed).\n'
    sys.stderr.write(msg_template.format(s=repr(seed)))
except:
    sys.stderr.write('An error interpreting the input occurred.  Expecting 3 arguments:\n')
    sys.stderr.write('    sample_size mean standard_deviation\n\nThe stacktrace of error:\n')
    raise

# This part is the actual simulation and output:
for i in xrange(n):
    variate = random.normalvariate(mu=mean, sigma=std_dev)
    print(variate)

