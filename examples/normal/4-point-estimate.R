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
}
