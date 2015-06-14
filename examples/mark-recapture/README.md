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
