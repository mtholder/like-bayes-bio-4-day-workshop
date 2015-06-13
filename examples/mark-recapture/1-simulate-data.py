#!/usr/bin/env python
'''Take at least 3 arguments:
    initial_n death_prob tag_loss_prob_1 [ tag_loss_prob_2  tag_loss_prob_3 ...]
writes the number of tags found in each year until 0 are seen
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
    death_prob = float(arguments[1])
    assert death_prob >= 0.0
    assert death_prob <= 1.0
    tag_loss_prob_list = [float(i) for i in arguments[2:]]
    for tlp in tag_loss_prob_list:
        assert tlp >= 0.0
        assert tlp <= 1.0
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
    sys.stderr.write('An error interpreting the input occurred.\n')
    sys.stderr.write(__doc__)
    raise

# This part is the actual simulation and output:
curr = n

while curr > 0:
    # set the tag loss prob to the next element (if there is one)
    if tag_loss_prob_list:
        tag_loss_prob = tag_loss_prob_list.pop(0)
    next_n = 0
    for i in xrange(curr):
        if random.random() < death_prob:
            continue
        if random.random() < tag_loss_prob:
            continue
        next_n = next_n + 1
    print(next_n)
    curr = next_n
        