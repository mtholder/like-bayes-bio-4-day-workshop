#!/usr/bin/env python

args <- commandArgs(trailingOnly=TRUE);
filename <- args[1]
n.iter <- as.integer(args[2]);
if (is.na(n.iter) || n.iter < 0) {
    write("The second argument should be the number of iterations", stderr());
    q(status=1);
}
initial.pop <- as.integer(args[3]);
if (is.na(initial.pop) || initial.pop < 0) {
    write("The third argument should be an initial population size", stderr());
    q(status=1);
}
tag.loss.prob.list <- as.numeric(args[4:length(args)]);
if (length(args) < 3) {
    write("The fourth argument must be tag loss probability");
    q(status=1);
}
## Read the data...
data <- as.vector(as.matrix(read.table(filename, header=FALSE)))

####################################################
# The following 3 functions vary from application to application...
####################################################
PROPOSAL.WINDOW.SIZE <- 0.01;
propose.new.state <- function(theta) {
    #Return a randomly proposed state and the  log(Hastings ratio).
    # propose a state from among the adjacent states
    u.neg.half.to.half <- runif(1) - .50;
    diff <- PROPOSAL.WINDOW.SIZE*u.neg.half.to.half;
    proposed.state <- theta + diff;
    if (proposed.state < 0.0) {
        proposed.state <- -proposed.state;
    }
    if (proposed.state > 1.0) {
        proposed.state <- 1.0 - (1.0 - proposed.state);
        proposed <- theta + 1;
    }
    return(c(proposed.state, 0.0));
}

calc.ln.likelihood <- function(data, death.prob) {
    ln.like <- 0.0;
    prev.n <- initial.pop;
    for (curr.index in 1:length(data)) {
        if (curr.index <= length(tag.loss.prob.list)) {
            tag.loss.prob <- tag.loss.prob.list[curr.index];
        }
        curr.n <- data[curr.index];
        n.exited <- prev.n - curr.n;
        recovery.prob <- (1.0 - tag.loss.prob)*(1 - death.prob);
        exit.prob <- 1.0 - recovery.prob;
        if (curr.n > 0) {
            ln.like <- ln.like + curr.n*log(recovery.prob);
        }
        if (n.exited > 0) {
            ln.like <- ln.like + n.exited*log(exit.prob);
        }
        prev.n <- curr.n;
    }
    return(ln.like);
}

calc.ln.prior = function(param) {
    return(0.0);
}

####################################################
# The code below is basically the same for any MCMC
####################################################
# This is MCMC using the Metropolis algorithm, except
# we often have to write out each visited state rather
# than keep a count of them in `mcmc.samples`
metropolis.hastings <- function(data, num.it, start.state, sample.freq) {
    state <- start.state;
    ln.likelihood <- calc.ln.likelihood(data, state);
    ln.prior <- calc.ln.prior(state);
    write("Iter\tlnL\tlnPrior\ttheta\n", stderr());
    mcmc.samples <- c();
    for (i in 1:num.it) {
        x <- propose.new.state(state);
        proposed.state <- x[1];
        ln.hastings.ratio <- x[2];
        proposed.ln.prior <- calc.ln.prior(proposed.state);
        proposed.ln.likelihood <- calc.ln.likelihood(data, proposed.state);
        
        ln.prior.ratio <- proposed.ln.prior - ln.prior;
        ln.likelihood.ratio <- proposed.ln.likelihood - ln.likelihood;

        ln.posterior.ratio <- ln.prior.ratio + ln.likelihood.ratio;
        ln.acceptance.ratio <- ln.posterior.ratio + ln.hastings.ratio;
        
        if (log(runif(1)) < ln.acceptance.ratio) {
            # Accept the move, update the state and current state info
            state <- proposed.state;
            ln.likelihood = proposed.ln.likelihood;
            ln.prior = proposed.ln.prior;
        }
        if ((i %% sample.freq) == 0) {
            # there are better ways in R than concatenating a vector
            #   in this way
            mcmc.samples <- c(mcmc.samples, state); 
            write(paste(i,
                        "\t",
                        ln.likelihood,
                        "\t",
                        ln.prior,
                        "\t",
                        state),
                 stderr());
        }
    }
    return(mcmc.samples)
}

start.state <- runif(1);
mcmc.samples <- metropolis.hastings(data, n.iter, start.state, 100);
mcmc.samples <- sort(mcmc.samples);
mean <- mean(mcmc.samples);
write(paste("posterior mean =", mean), stdout());
lower.ci.index <- as.integer(0.025*length(mcmc.samples));
upper.ci.index <- as.integer(0.975*length(mcmc.samples));
lower.ci <- mcmc.samples[lower.ci.index];
upper.ci <- mcmc.samples[upper.ci.index];

write(paste("95% credible interval [", 
            lower.ci,
            ", ",
            upper.ci,
            "]"),
      stdout());

