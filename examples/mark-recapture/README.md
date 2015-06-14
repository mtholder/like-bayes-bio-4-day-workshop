# simplified mark recapture problem

Can we estimate the per-month probability of death from 
tag + recapture. 

Here we assume that the tags have transmitters, so that the
animals can be found.

Some number N are tagged at the beginning of the study.

Every month you survey the area. Assume that you can detect 
every functioning tag attached to a living organism.

You know (from other studies) that the probabilities that a tag
stops functioning (or falls off) are:
  * 0.10 in the first month,
  * 0.15 in the second month,
  * 0.2 in the third month, and
  * 0.25 for every month after that.

We want to know the per month death probability. We assume that this
number is constant across time and across individuals.

The tricky thing is that both tag failure and death contribute to 
dropping counts.  And the tag loss probabilities are quirky. 

# data set simulation

    python 1-simulate-data.py 1000 .1 .1 .15 .2 .25 >data.tsv 2>seed-log.txt

created the data. 1000 is the number of individuals tagged.
The first .1 is the per-month death probability.
The remaining parameters are the per-month tag loss probabilities.


# estimation:
First run

    python do-it-yourself-interval-estimation.py data.tsv 1000 .1 .15 .2 .25

It will crash becuase it does not have the per-datum log likelihood calculation
filled in. See if you can figure out what the log-likelihood formula is.

*hint* treat the fate of the tagged individuals from the previous month 
  (the fate can be either "detected again" or "exited") as independent for 
  each month.


If you work that formula for the log-likelihood, see if you can program it on
line 36 of the code.

The comments between "DO NOT EDIT ABOVE HERE" and "DO NOT EDIT BELOW HERE"
explain how the variables that are in scope at this spot in the code correspond
the the quantities that are relevant to the likelhood calculation.


# Seeing a worked example of interval estimation:
If you prefer, you can see how I solved the problem using:

    python 5-interval-estimate.py data.tsv 1000 .1 .15 .2 .25

As you would expect to happen 95% of the time, the 95% confidence interval contains
the true "death probability" that I used in the simulations (it was 0.1)
