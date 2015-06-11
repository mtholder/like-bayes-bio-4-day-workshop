#!/usr/bin/env R
# Take 3 arguments:
#    sample_size mean standard_deviation
# and writes to standard output a simulated realization from the normal distribution
# (one variate per line).

# If you were trying to do this in one line, you could use:
#    write(paste(rnorm(n, mean, std.dev), sep="\n"), stdout());

# In a command-line program, we have 2 output "streams" that we can write to
# You should used stderr for status messages and error messages; and
#   use stdout for the "real" output of the program. "write" writes to 
#   stdout or stderr() depending on the argument.
# Some notes on the output below:
#   Writing to the ouput stream requires that we convert our value to a "string"
#   We  can join different fragments to be written using paste with a c()
#       construct to describe the series of values. The default separator is 
#       a space.
#   The "\n" represents a newline character.
args <- commandArgs(trailingOnly=TRUE);
n <- as.integer(args[1]);
if (is.na(n) || n < 1) {
    write("The first argument should be a sample size", stderr());
    q(status=1);
}
mean <- as.numeric(args[2]);
if (is.na(mean)) {
    write("The second argument should be a mean", stderr());
    q(status=1);
}
std.dev <- as.numeric(args[3]);
if (is.na(std.dev) || std.dev < 0) {
    write("The third argument should be a standard deviation", stderr());
    q(status=1);
}

# Whenever you write a program with a random number usage, it is wise to:
#   1. report the pseudorandom number generator's seed, and
#   2. supply some way of specifying the seed.
# Here we use the "environment" that the script is run in. This is cryptic
#   (not listed as a command line argument), but unobtrusive.  And we don't
#   expect to use this option often in this script
seed <- Sys.getenv('SIMULATE_DATA_SET_SEED');
if (seed == "") {
    seed <- as.numeric(Sys.time());
} else {
    seed <- as.numeric(seed);
    if (is.na(seed)) {
        write("The SIMULATE_DATA_SET_SEED value should be an integer", stderr());
        q(status=1);
    }
}
set.seed(seed);
write(paste("pseudorandom number seed =",
              seed,
              "(use the env. var. SIMULATE_DATA_SET_SEED to set the seed)."),
      stderr());

# This part is the actual simulation and output:
realization <- rnorm(n, mean, std.dev);
write(paste(realization, sep="\n"), stdout());
