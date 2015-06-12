#!/usr/bin/env python
# This program 2 or 3 arguments:
#      a filename of a file with the data (one datum per line), and
#      a standard deviation
# Reports summary statistics and an MLE for mu
# If a third argument is provided, it is interpreted as a value for mu
#     and the likelihood is calculated.

args <- commandArgs(trailingOnly=TRUE);
filename <- args[1]
std.dev <- as.numeric(args[2]);
if (is.na(std.dev) || std.dev < 0) {
    write("The second argument should be a standard deviation", stderr());
    q(status=1);
}
user.mu <- NA;
if (length(args) > 2) {
    user.mu <- as.numeric(args[3]);
    if (is.na(user.mu)) {
        write("The third argument must be a value for the mean");
        q(status=1);
    }
}
## Read the data...
data <- as.vector(as.matrix(read.table(filename, header=FALSE)))
avg = mean(data);
print(paste("sample mean =", avg));

calc.ln.likelihood = function(data, mu) {
    coefficient <- 1.0/(std.dev * sqrt(2*pi));
    ln.c <- log(coefficient);
    var.factor <- -1.0/(2.0*std.dev*std.dev);
    diff.squared <- (data - mu)^2;
    datum.ln.like <- ln.c + var.factor*(data - mu)^2;
    return(sum(datum.ln.like));
}

optimizer.ln.like = function(param) {
    return(-calc.ln.likelihood(data, param));
}

find.lower.bound <- function(data,
                             ln.L.func,
                             mle,
                             ln.L.cutoff) {
    eps <- 0.1; # to avoid crashes if the MLE= 0
    lb <- mle*.95 - eps;
    lb.ln.L <- calc.ln.likelihood(data, lb);
    while (lb.ln.L > ln.L.cutoff) {
        value.diff <- mle - lb;
        lb <- lb - value.diff;
        lb.ln.L <- ln.L.func(data, lb);
    }
    return(lb);
}
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
    lb <- find.lower.bound(data, ln.L.func, mle, ln.L.cutoff)
    return(optimize(f = ln.likelihood.lower.conf,
                    interval=c(lb, mle))$minimum);
}

find.upper.bound <- function(data,
                             ln.L.func,
                             mle,
                             ln.L.cutoff) {
    eps <- 0.1; # to avoid crashes if the MLE= 0
    ub <- mle*1.05 + eps;
    ub.ln.L <- calc.ln.likelihood(data, ub);
    while (ub.ln.L > ln.L.cutoff) {
        value.diff <- ub - mle;
        ub <- ub + value.diff;
        ub.ln.L <- ln.L.func(data, ub);
    }
    return(ub);
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
    ub <- find.upper.bound(data, ln.L.func, mle, ln.L.cutoff)
    return(optimize(f = ln.likelihood.upper.conf,
                    interval=c(ub, mle))$minimum);
}

if (!is.na(user.mu)) {
    user.ln.like <- calc.ln.likelihood(data, user.mu);
    print(paste("ln[Pr(data | mu={m})] = ", user.ln.like));
} else {
    solution <- optimize(f=optimizer.ln.like,
                         interval=c(min(data), max(data)));
    mle <- solution$minimum;
    max.ln.like <- -solution$objective;
    print(paste("ln[Pr(data | mu=",
                mle,
                ")] =",
                max.ln.like));
    chi.sq.critical <- 3.84
    cutoff.diff <- chi.sq.critical/2.0
    cutoff <- max.ln.like - cutoff.diff
    lower.mu = find.lower.confidence.bound(data,
                                           ln.L.func=calc.ln.likelihood,
                                           mle=mle,
                                           ln.L.cutoff=cutoff)
    lower.ln.l = calc.ln.likelihood(data, lower.mu)
    print(paste("lower-CI ln[Pr(data | lower.mu=",
                 lower.mu,
                 "})] = ",
                 lower.ln.l));
    upper.mu = find.upper.confidence.bound(data,
                                           ln.L.func=calc.ln.likelihood,
                                           mle=mle,
                                           ln.L.cutoff=cutoff)
    upper.ln.l = calc.ln.likelihood(data, upper.mu)
    print(paste("upper-CI ln[Pr(data | lower.mu=",
                 upper.mu,
                 "})] = ",
                 upper.ln.l));
    print(paste("Approx. 95% CI: (",
                 lower.mu,
                 ",",
                 upper.mu,
                 ")"));
}
