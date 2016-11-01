# Round to integer keeping the sum fixed
# Code taken from roundfixS function of sfsmisc-package


#' Round to integer keeping the sum fixed
#'
#' @param x A numeric vector which must sum to an integer
#'
#' @return A vector of rounded integers summing up to sum(x).
#' @references Code taken from roundfixS function of sfsmisc-package by Martin
#'   Maechler
#' @export
round_fix_sum <- function(x) {
   n <- length(x);
   x0 <- floor(x);
   e <- x - x0;
   S. <- sum(e);
   stopifnot(all.equal(S., (S <- round(S.))));
   if (S > 0) {
      r <- numeric(n);
      r[sort.list(e, decreasing = TRUE)[1:S]] <- 1;
      x0 <- x0 + r;
      return(x0)
   } else {
      return(x);
   }
}
