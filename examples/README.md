# code examples for likelihood + Bayesian approaches course

Examples are in `python` and `R`. The python demos depend on `scipy`. See
http://mtholder.github.io/like-bayes-bio/
for a discussion of installing these dependencies.

  1. the `normal` subdirectory demonstrates how a simple problem of estimation
    and/or testing an *a priori* value for a parameter can be done using ML.

  2. `mark-recapture` shows how the exact same techniques can be used in a
   slightly more realistic (but still remarkably simplistic) example. 
   What is cool is that the scripts for interval estimation in this 
   directory and the `normal` dir only differ in the calculation of the
   likelihood, and a slight tweak to make the optmization code obey the 
   bounds on a parameter that is a probability. Everything else is reusable.