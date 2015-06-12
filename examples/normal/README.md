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

We can pass in the simulated data to `2-calc-likelihood.py` and also give that script a
standard deviation and a guess for the population mean:

    python 2-calc-likelihood.py data.tsv 1 10

would report the likelihood at the true (simulated) value. The script also reports
the sample mean for the data set.

See if you can find a better estimate of the population mean by calling the script 
several times and paying attention to the reported log-likelihood.

## numerical stability

Note that `2-calc-likelihood.py` calculates a likelihood and then simply takes the log
to get the log likelihood. Why is this a bad idea?

If you try to score a very implausible parameter value:

    python 2-calc-likelihood.py data.tsv 1 100 

you'll run into problems because the likelihood is so small that it can no longer
be represented by a naive computer program without running into problems associated
with [floating point underflow](https://en.wikipedia.org/wiki/Arithmetic_underflow).

Note how the code is slightly different in `3-stable-calc-likelihood.py`, and this 
allows us to get a likelihood:
    
    python 3-stable-calc-likelihood.R data.tsv 1 100 

or even much more extreme values.


## Point estimation
Calculating the likelihood is nice for testing, but usually we more interested in
estimation.

     python 4-point-estimate.py data.tsv 1

will estimate a point estimate - the MLE of the mean.  We'll discuss the computational
methods used by the script a little bit in the second lecture. Basically, we pass
a function to the `optimize.brent` (in Python) or `optimize` (in R) functions.
They try out lots of parameter values (in an intelligent manner), and return the 
parameter value that gave the lowest score.
Since we want to **maximize** the likelihood, we insert a tiny function that takes
or log-likelihood and turns it into a negative log-likelihood. That is the function
that we hand to the optimizer. The minimum of the negative log likelihood will occur
at the parameter value that produces the maximum likelihood.


It you pass in a parameter value for the mean:

     python 4-point-estimate.py data.tsv 1 10

the script will use a likelihood ratio test to tell you if you can reject that value of the
parameter.

## Interval estimation
We can use the chi-squared distribution's cutoff to produce a 95% confidence interval.

    python 5-interval-estimate.py data.tsv 1

will report this interval estimate (as well as the interval that you would get if you
recognized that this inference problem was a perfect case for using a 
[**Z**-test](https://en.wikipedia.org/wiki/Z-test) instead of using ML methods).

Note that we just use the same optimization routines to find the boundary points, but we
just tweak the function to be optimized. We want to find the parameter where the lnL is
3.84/2 lower than the maximum lnL. So we just use the absolute value of the difference between
the lnL and the maxLnL as our objective function.

We also have to used bounded optimization methods (because the lnL will cross that threshold 
value at least twice: once below the MLE and once above).



   
