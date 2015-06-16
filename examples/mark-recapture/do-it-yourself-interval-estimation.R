#!/usr/bin/env python
# This program 2 or 3 arguments:
#      a filename of a file with the data (one datum per line), and
#      a standard deviation
# Reports summary statistics and an MLE for prob
# If a third argument is provided, it is interpreted as a value for prob
#     and the likelihood is calculated.

args <- commandArgs(trailingOnly=TRUE);
filename <- args[1]
initial.pop <- as.integer(args[2]);
if (is.na(initial.pop) || initial.pop < 0) {
    write("The second argument should be an initial population size", stderr());
    q(status=1);
}
tag.loss.prob.list <- as.numeric(args[3:length(args)]);
if (length(args) < 3) {
    write("The third argument must be tag loss probability");
    q(status=1);
}
## Read the data...
data <- as.vector(as.matrix(read.table(filename, header=FALSE)))
avg = mean(data);

calc.ln.likelihood = function(data, death.prob) {
    ln.like <- 0.0;
    prev.n <- initial.pop;
    for (curr.index in 1:length(data)) {
        if (curr.index <= length(tag.loss.prob.list)) {
            tag.loss.prob <- tag.loss.prob.list[curr.index];
        }
        curr.n <- data[curr.index];
        n.exited <- prev.n - curr.n;
        ###########################################
        # DO NOT EDIT ABOVE HERE
        # curr.n is the number of tags seen in this month.
        # n.exited is the difference between this number and the previous month.
        # tag.loss.prob is the probability of a tag that is working at the 
        #   beginning of the previous month falls off or stops working during
        #   the previous month
        # death.prob is the probability that an individual alive at the 
        #   beginning of the previous month is dead at this point in time.
        ###########################################
        if (curr.n > 0) {
            ln.like <- ln.like + INSERT CODE HERE
        }
        if (n.exited > 0) {
            ln.like <- ln.like + INSERT CODE HERE
        }
        ###########################################
        # DO NOT EDIT BELOW HERE
        ###########################################
        prev.n <- curr.n;
    }
    q(status=1)
    return(ln.like);
}

optimizer.ln.like = function(param) {
    return(-calc.ln.likelihood(data, param));
}

EPSILON <- 1e-6;
find.lower.confidence.bound <- function(data,
                                ln.L.func,
                                mle,
                                ln.L.cutoff) {
    ln.likelihood.lower.conf <- function(x) {
        if (x > mle) {
            return(Inf);
        }
        ln.L <- calc.ln.likelihood(data, x)
        return(abs(ln.L - ln.L.cutoff));
    }
    return(optimize(f = ln.likelihood.lower.conf,
                    interval=c(0, mle))$minimum);
}

find.upper.confidence.bound <- function(data,
                                ln.L.func,
                                mle,
                                ln.L.cutoff) {
    ln.likelihood.upper.conf <- function(x) {
        if (x < mle) {
            return(Inf);
        }
        ln.L <- calc.ln.likelihood(data, x)
        return(abs(ln.L - ln.L.cutoff));
    }
    return(optimize(f = ln.likelihood.upper.conf,
                    interval=c(mle, 1.0 - EPSILON))$minimum);
}

solution <- optimize(f=optimizer.ln.like,
                     interval=c(0, 1 - EPSILON));
mle <- solution$minimum;
max.ln.like <- -solution$objective;
print(paste("ln[Pr(data | prob=",
            mle,
            ")] =",
            max.ln.like));
chi.sq.critical <- 3.84
cutoff.diff <- chi.sq.critical/2.0
cutoff <- max.ln.like - cutoff.diff
lower.prob = find.lower.confidence.bound(data,
                                       ln.L.func=calc.ln.likelihood,
                                       mle=mle,
                                       ln.L.cutoff=cutoff)
lower.ln.l = calc.ln.likelihood(data, lower.prob)
print(paste("lower-CI ln[Pr(data | lower.prob=",
             lower.prob,
             "})] = ",
             lower.ln.l));
upper.prob = find.upper.confidence.bound(data,
                                       ln.L.func=calc.ln.likelihood,
                                       mle=mle,
                                       ln.L.cutoff=cutoff)
upper.ln.l = calc.ln.likelihood(data, upper.prob)
print(paste("upper-CI ln[Pr(data | lower.prob=",
             upper.prob,
             "})] = ",
             upper.ln.l));
print(paste("Approx. 95% CI: (",
             lower.prob,
             ",",
             upper.prob,
             ")"));
