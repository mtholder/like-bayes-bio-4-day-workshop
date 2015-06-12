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
print(paste("mean =", avg));

calc.likelihood = function(data, mu) {
    datum_likelihood_vec <- dnorm(data, mean=mu, sd=std.dev);
    return(prod(datum_likelihood_vec));
}

naive.calc.ln.likelihood = function(data, mu) {
    return(log(calc.likelihood(data, mu)));
}

calc.ln.likelihood = naive.calc.ln.likelihood;
if (!is.na(user.mu)) {
    user.ln.like <- calc.ln.likelihood(data, user.mu);
    print(paste("ln[Pr(data | mu={m})] = ", user.ln.like));
}
