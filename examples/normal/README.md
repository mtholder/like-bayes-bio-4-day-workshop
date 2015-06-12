# Toy example of estimating the mean of a normal distribution.

## simulating a data set

    python 1-simulate-data.py 20 10 1

or 

    Rscript 1-simulate-data.R 20 10 1

will print out a simple dataset of 20 realizations from a normal distribution
with a mean of 10.0 and a standard deviation of 1.0.
Both scripts demonstrate strategies for error reporting, so you can try to
pass in nonsensical values for the arguments and see if the errors are reported
in a reasonable way.

The scripts also have some discussion of
  1. making sure that your random-number-dependent scripts give the user
   a means of setting a seed for the pseudorandom number generator; and
  2. using both the standard "error" and standard "output" stream for communicating
   to the user.

I generated the data with:

    python 1-simulate-data.py 20 10 1 >data.tsv 2>sim-log.txt

to create the simulated data in `data.tsv` and the status messages (including the seed
used in `sim-log.txt`).


## calculating a likelihood.

