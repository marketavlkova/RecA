#!/usr/local/bin/Rscript
### function to find all divisors

divisors <- function(x) {
  y <- seq_len(x)
  y[x%%y == 0]
}
