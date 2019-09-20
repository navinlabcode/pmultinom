## Functions used internally
# These functions don't check their input; it's up to the user-facing functions
# to make sure input is sensible

## Reference implementation for testing
reference.cdf <- function(upper, probs, size, method=NULL)
{
  if (length(upper)==0) return(1)
  stopifnot(all.equal(sum(probs), 1))
  # Probability of filling the first cell
  first.cell.dist <- stats::dbinom(0:min(upper[1], size), size, probs[1])
  # Probability of filling the remaining cells
  conditional.cdf <- sapply(0:min(upper[1], size), function(m) reference.cdf(upper[-1], probs[-1]/sum(probs[-1]), size-m))
  return(sum(first.cell.dist * conditional.cdf))
}

reference.complementary.cdf <- function(lower, probs, size, method=NULL)
{
  if (length(lower)==0) return(1)
  stopifnot(all.equal(sum(probs), 1))
  if (lower[1]+1 > size) return(0)
  # Probability of filling the first cell
  first.cell.dist <- stats::dbinom((lower[1]+1):size, size, probs[1])
  # Probability of filling the remaining cells
  conditional.cdf <- sapply((lower[1]+1):size, function(m) reference.complementary.cdf(lower[-1], probs[-1]/sum(probs[-1]), size-m))
  return(sum(first.cell.dist * conditional.cdf))
}

# Truncated Poisson pmf
# Seems straightforward, but turns out to be tricky to implement while
# maintaining numerical stability when the event we're conditioning on has very
# small probability
tpois.pmf <- function(xs, a, b, l)
{
  # If log(l) is -Inf, we'll use the following reasoning to fill in the answer
  # with a limit. If a Poisson random variable has to be above a, but it is so
  # small that it is extremely unlikely that it exceeds this condition, then it
  # is even more unlikely that it exceeds it by more than the bare minimum. So
  # it has to be a+1.
  # This is also the case if b == a + 1
  if ((identical(log(l), -Inf) & a >= 0) | (b == a+1)) {
    return(ifelse(xs == a + 1, 1, 0))
  } else {
    log.num <- stats::dpois(xs, l, log=TRUE)

    # If l is 0, but a >= 0, the condition will have probability 0. For l values
    # very close to zero, I use the following formula, which is accurage for
    # small l. Actually it doesn't have to be that close to zero.
    # I should get an error bound and use that in the fuutre.
    log.denom <- log(stats::ppois(b, l) - stats::ppois(a, l))
    if (is.infinite(log.denom)) log.denom <- (a+1) * log(l) - log(a+1) - lgamma(a+1)

    return(ifelse(a < xs & xs <= b,
                  exp(log.num - log.denom),
                  0))
  }
}

# Functions to solve P(X1>b1, X2>b2, ... | Lambda1=lambda*p1,
# Lambda2=lambda*p2, ...) in the Poisson approximation. First, for inverting
# when there's only one condition, then, for inverting when there's multiple
# conditions.
#' @useDynLib pmultinom sequentialNewtonRootR
invert.equiprobable <- function(n.categories, target.probability, cutoff)
{
  # Solve P[X>=cutoff]^n.categories = target.probability, assuming a Poisson
  # distribution, solving for the Poisson rate
  qgamma(target.probability^(1/n.categories), shape = cutoff, rate = 1)
}
invert.poisson.survival.function <- function(ps, bs, ms, target, eps=10^-2)
{
  stopifnot(length(ps) == length(bs) & length(ps) == length(ms))
  
  # bs are the highest excluded value
  # The lowest included value:
  cs <- bs + 1
  
  # Figure out how hard each condition is to satisfy, quantified as how many
  # individuals it would require
  indiv.lambdas <- invert.equiprobable(ms, target, cs) / ps
  
  # If there's only one condition, this can be solved exactly using R's builtin
  # functionality for inverting the incomplete gamma function
  if (length(indiv.lambdas) == 1)
  {
    return(indiv.lambdas)
  }
  
  # Sort from hardest to easiest to satisfy
  ordering <- order(indiv.lambdas, decreasing=TRUE)
  ps.sorted <- ps[ordering]
  cs.sorted <- cs[ordering]
  ms.sorted <- ms[ordering]
  
  # Find lambda that achieves the target
  # C call for sequentially finding Poisson rates to satisfy prefixes of the conditions:
  ccall <- .C("sequentialNewtonRootR", PACKAGE="pmultinom",
    npop=length(indiv.lambdas), ps=as.double(ps.sorted),
    cs=as.integer(cs.sorted), ms = as.integer(ms.sorted),
    lambdaStart = as.double(indiv.lambdas[ordering[1]]),
    target=as.double(target), eps=as.double(eps), output=double(1))
  # And the final result:
  return(ccall$output)
}
# Functions for the Fourier convolution method.  I shortsightedly gave these
# names like "sum.tpois.pmf" as if there was only one way to calculate that.

# Convolve function
# Inspired by "convolve" in the stats package
fftw.convolve <- function(unpadded.x, unpadded.y, plan)
{
  # Assuming without checking that length(x) == length(y)
  n.in <- length(unpadded.x)
  n.padded <- 2*length(unpadded.x) - 1
  x <- c(rep.int(0, n.in - 1), unpadded.x)
  y <- c(unpadded.y, rep.int(0, n.in - 1))
  complex.result <- fftw::IFFT(fftw::FFT(x,plan=plan) * Conj(fftw::FFT(y,plan=plan)), plan=plan)
  # Assume without checking that the input was real and that the output,
  # therefore, should also be real
  Re(complex.result)
}

# Probability mass function of the sum of truncated Poisson distributions.
# Calculated in a pretty naive and general way: calculate the pmfs of the
# individual truncated Poisson random variables and convolve them.
sum.tpois.pmf <- function(as, bs, ls, n)
{
  # "bs" are max excluded, ls are proportional to clone frequencies
  # Calculate a matrix, where each column is a single truncated Poisson pmf
  xs <- 0:n
  tpois.pmfs <- mapply(function(a, b, l)
    tpois.pmf(xs, a, b, l),
    as, bs, ls, SIMPLIFY=TRUE)
  # Take the convolutions. Note that the convolve function needs the second
  # argument to be reversed for some reason (see the help page of
  # stats::convolve), which is accomplised here by flipping the rows of the
  # matrix
  plan <- fftw::planFFT(nrow(tpois.pmfs)*2-1, effort=0)
  sum.tpois.pmf <- Reduce(function(x, y) fftw.convolve(x, y,plan=plan)[(0:n)+1],
                          data.frame(tpois.pmfs[nrow(tpois.pmfs):1,-1]),
                          tpois.pmfs[,1])
  return(sum.tpois.pmf)
}

# Calculates the multinomial probability, given the distribution of the
# appropriate sum of truncated Poisson random variables
prob.between <- function(as, bs, ps, n, sum.pmf, l)
{
  ls <- l * ps
  # Set negative values to zero before logging. In practice these are only very
  # slightly negative so its fine.
  log.cond.prob.n <- log(pmax(sum.pmf[n + 1], 0))
  log.uncond.prob.n <- stats::dpois(n, sum(ls), log=TRUE)
  log.poisson.event.prob <- sum(log(stats::ppois(bs, ls) - stats::ppois(as, ls)))
  return(exp(log.poisson.event.prob + log.cond.prob.n - log.uncond.prob.n))
}

# Functions for the Edgeworth method
# Calculate a cumulant of the truncated Poisson distribution
#' @useDynLib pmultinom calculate_correction
# This can be like half the computational cost, only gives eleven cumulants,
# and only seven or so of them accurately, so I need to replace this
cumulant <- Vectorize(function(lambda, cutoff, cumulant)
{
  # Special case for the mean. I'd fill in lambda=0 special cases for all of
  # them if I knew them
  if (lambda == 0 & cumulant == 1) return(cutoff)
  stopifnot(cumulant <= 11)
  x <- exp(-lambda + cutoff*log(lambda) - (lgamma(cutoff) + pgamma(lambda, shape=cutoff, rate=1, log=TRUE)))
  correction <- .C("calculate_correction", PACKAGE="pmultinom",
     cumulant=as.integer(cumulant), c=as.integer(cutoff),
     ld=as.double(lambda), x = as.double(x), correction=double(1))$correction
  return(lambda + correction)
})
# Calculate probabilist's Hermite numbers
# This is a small but significant chunk of the computation time in the
# Edgeworth approximation calculations, and it can easily be calculated much
# faster recursively, so that's a to-do
prob.hermite.number <- function(n)
{
  Re((1i)^n * (1/2)*gamma((1/2)*n + 1/2)*2^((1/2)*n)*(1 + (-1)^n)/sqrt(pi))
}
# "advance" function from PDQutils
# Slows it down a lot, so I'll have to find another way to do this
advance <- function (kms) 
{
    current <- kms$current
    mold <- kms$mold
    n <- length(current)
    ords <- 1:(length(current))
    stopifnot(n == sum(current * ords))
    sumcur <- n
    m <- 1
    is.done <- FALSE
    while (!is.done) {
        sumcur <- sumcur - current[m] * m + (m + 1)
        current[m] <- 0
        current[m + 1] <- current[m + 1] + 1
        m <- m + 1
        is.done <- (sumcur <= n) || (m > mold)
    }
    mold <- max(mold, m)
    current[1] <- n - sumcur
    retv <- list(current = current, mold = mold)
    return(retv)
}
# PMF for a sum of truncated Poisson distributions, using an Edgeworth
# expansion. A small but significant fraction of the computation time is
# within this function, so a lot of this looping should be moved to the C
# code. Also, this function is misleadingly named. It is not a probability
# mass function, and in fact only provides the probability mass at the mean.
tpois.sum.pmf.edgeworth <- function(lambdas, cutoffs, multiplicities, eps=10^-5)
{
  combined.mean <- sum(multiplicities * cumulant(lambdas, cutoffs, 1))
  combined.variance <- sum(multiplicities * cumulant(lambdas, cutoffs, 2))
  cumulants <- c(0, 1)
  already.calculated <- 2
  current <- 3
  Sn <- numeric(0)
  coef_vec <- numeric(0)
  poly_vec <- integer(0)
  contribution <- 0
  last.contribution <- 0
  for (s in 1:1000) {
    # Calculate the current cumulant
    current.cumulant.index <- s + 2
    Sn <- c(Sn,
        sum(multiplicities * cumulant(lambdas, cutoffs, current.cumulant.index)) /
        sqrt(combined.variance)^current.cumulant.index
    )
    nexterm <- 0
    kms <- list(mold = 1, current = c(s, rep(0, s - 1)))
    last.contribution <- contribution
    contribution <- 0
    while (kms$mold <= s) {
      r <- sum(kms$current)
      coefs <- ((Sn[1:s]/factorial(3:(s + 2)))^kms$current)/factorial(kms$current)
      coef <- prod(coefs)
      # Record the coefficient, and the index of the Hermite polynomial it's a
      # coefficient on
      coef_vec <- c(coef_vec, coef)
      poly_vec <- c(poly_vec, s + 2 * r)
      # Add to the contribution to decide whether to stop
      contribution <- contribution + coef * prob.hermite.number(s + 2 * r)
      kms <- advance(kms)
    }
    # Decide whether to stop
    if (s > 4 & abs(contribution) < eps & abs(last.contribution) < eps) break
    if (s >= 1000) stop("Edgeworth expansion did not converge")
  }
  # Calculate the correction
  correction <- sum(coef_vec * prob.hermite.number(poly_vec))
  # Calculate the probability
  (1 / sqrt(combined.variance)) * dnorm(0) * (1 + correction)
}

## Exported functions

#' Calculate the probability that a multnomial random vector is between,
#' elementwise, two other vectors.
#'
#' @param lower Vector of lower bounds. Lower bounds are excluded
#' @param upper Vector of upper bounds. Upper bounds are included
#' @param size Number of draws
#' @param probs Cell probabilities
#' @param method Method used for computation. Only method currently implemented is "exact"
#' @return The probability that a multinomial random vector is greater than all
#'   the lower bounds, and less than or equal all the upper bounds:
#'
#' P(l1 < N1 <= u1, ..., lk < Nk <= uk)
#'
#' If only the upper bounds are given, then this is the multinomial CDF:
#'
#' P(N1<=u1, ..., Nk<=uk)
#'
#' If only the upper bounds are given, then this is the multinomial tail probability:
#'
#' P(N1>l1, ..., Nk>lk)
#'
#' @details
#' The calculation follows the scheme suggested in Levin (1981): begin with the
#' equivalent probability for a Poisson random vector, and update it by
#' conditioning on the sum of the vector being equal to the size parameter,
#' using Bayes' theorem. This requires computation of the distribution of a sum
#' of truncated Poisson random variables, which is accomplished using
#' convolution, as per Levin's suggestion for an exact calculation. Levin's
#' suggestion for an approximate calculation, using Edgeworth expansions, may be
#' added to a later version. Fast convolution is achieved using the fastest
#' Fourier transform in the west (Frigo, Johnson 1998).
#' 
#' @examples
#' # To determine the bias of a die, Rudolph Wolf rolled it
#' # 20,000 times. Side 2 was the most frequently observed, and
#' # was observed 3631 times. What is the probability that a
#' # fair die would have a side observed this many times or
#' # more?
#'
#' # Input: 
#' 1 - pmultinom(upper=rep.int(3630, 6), size=20000,
#'               probs=rep.int(1/6, 6), method="exact")
#' # Output:
#' # [1] 7.379909e-08
#' 
#' # Therefore we conclude that the die is biased. Fougere
#' # (1988) attempted to account for these biases by assuming
#' # certain manufacturing errors. Repeating the calculation
#' # with the distribution Fougere derived:
#' 
#' # Input:
#' theoretical.dist <- c(.17649, .17542, .15276, .15184, .17227, .17122)
#' 1 - pmultinom(upper=rep.int(3630, 6), size=20000,
#'               probs=theoretical.dist, method="exact")
#' # Output:
#' # [1] 0.043362
#' 
#' # Therefore we conclude that the die still seems more biased
#' # than Fougere's model can explain.
#'
#' @references
#'   Fougere, P. F. (1988). Maximum entropy calculations on a discrete probability space. In Maximum-Entropy and Bayesian Methods in Science and Engineering (pp. 205-234). Springer, Dordrecht. doi:10.1007/978-94-009-3049-0_10
#' 
#'   Frigo, Matteo, and Steven G. Johnson. (2005). The design and implementation of FFTW3. Proceedings of the IEEE, 93(2), 216-231. doi:10.1109/JPROC.2004.840301
#'
#'   Levin, Bruce. (1981). "A Representation for Multinomial Cumulative Distribution Functions". Annals of Statistics 9 (5): 1123–6. doi:10.1214/aos/1176345593
#' 
#' @seealso \code{\link{invert.pmultinom}}
#'
#' @export
pmultinom <- function(lower=-Inf, upper=Inf, size, probs, method)
{
  if (!isTRUE(method %in% c("exact", "edgeworth"))) stop('Method must be specified. Only available method right now is "exact"')
  # Checking input
  stopifnot(all(is.numeric(lower) | is.na(lower)))
  stopifnot(all(is.numeric(upper) | is.na(upper)))
  stopifnot(all(is.numeric(size) | is.na(size)))
  stopifnot(all(is.numeric(probs) | is.na(probs)))
  stopifnot(all(sapply(0 <= probs & probs <= 1, isTRUE) | is.na(probs)))
  stopifnot(all(is.na(size)) | isTRUE(all.equal(size[!is.na(size)], round(size[!is.na(size)]))))
  # Test that the arguments recycle properly
  tryCatch(mapply(function(x, y, z) NULL, lower, upper, probs),
           warning=function(w) stop(sprintf("Error generated from warning in %s: %s",
                                            format(conditionCall(w)), conditionMessage(w))))
  # If the all of the sizes are missing, return NA
  if (all(is.na(size))) return(NA)
  # Make recycled versions of the arguments
  r.lower <- mapply(function(x, y, z) x, lower, upper, probs)
  r.upper <- mapply(function(x, y, z) y, lower, upper, probs)
  r.probs <- mapply(function(x, y, z) z, lower, upper, probs)
  # As much as I hate to add a special case for this, I need to check for upper
  # > lower. The problem is that the conditional distribution this algorithm
  # relies on is not defined otherwise. I do this check before some of the NA
  # checks because mathematically impossible conditions actually make a
  # probability possible to calculate even when inputs are missing.
  if (any(sapply(r.lower >= r.upper, isTRUE)))
  {
    return(0)
  }
  # If any of the probabilities are missing, return NA
  # I do this before checking normalization, because if some are missing, they
  # will not sum to 1, but the correct answer is still NA
  if (any(is.na(probs))) return(NA)
  # Check that the probabilities sum to 1, after recycling
  if (!isTRUE(all.equal(sum(r.probs), 1))) stop(sprintf("Probabilities do not sum to 1 after recycling: %s", paste(r.probs, collapse=", ")))
  # If any of the bounds are missing, return NA
  # I do this after checking normalization, because if the probabilities don't
  # normalize, the output is nonsense, no matter what the bounds are
  if (any(is.na(c(lower, upper)))) return(NA)

  # Convolution-based method
  if (method=="exact")
  {
    # Current method is to calculate all values from 0 to the largest input size.
    # An easy optimization would be to skip the ones below the smallest input
    # size.
    # The algorithm has a tuning parameter, and for each value of the tuning
    # parameter, is numerically stable for a range of of input values. Therefore,
    # the following loop calculates different parts of the path 0:n with different
    # values of the tuning parameter.
    a <- 3
    v <- a
    n <- v^2 - 1
    distribution <- sum.tpois.pmf(r.lower, r.upper, n * r.probs, n)
    path <- prob.between(r.lower, r.upper, r.probs, 0:n, distribution, n)
    while (n < max(size, na.rm=TRUE))
    {
      old.n <- n
      v <- v + a
      n <- v^2 -1
      distribution <- sum.tpois.pmf(r.lower, r.upper, n * r.probs, n)
      path <- c(path, prob.between(r.lower, r.upper, r.probs, (old.n+1):n, distribution, n))
    }
    return(path[size+1])
  } else if (method == "edgeworth")
  # Method based on an Edgeworth expansion
  {
    if (!all(upper == Inf)) stop("Edgeworth expansion only implemented for lower bounds so far")
    # Need to choose lambda to center the mean of truncated poissons on
    # "size", to make sure an edgeworth expansion will be good. This requires
    # solving an equation, and I'm choosing an arbitrary range to search for a
    # solution in, for now. This needs to be fixed to make it robust
    poisson.rate <- uniroot(function(x) sum(cumulant(x*probs, lower+1, 1)) - size, lower=0, upper=2*size)$root
    exp(
        # Prior
        sum(ppois(lower, poisson.rate*probs, lower.tail=FALSE, log=TRUE)) +
        # Probability of data under hypothesis
        log(tpois.sum.pmf.edgeworth(poisson.rate*probs,
            lower+1, rep.int(1,length(probs)))) -
        # Overall probability of data
        dpois(size, poisson.rate, log=TRUE)
    )
  }
}

#' Calculate the sample size such that the probability of a result is a given amount.
#'
#' @param lower Vector of lower bounds. Lower bounds are excluded
#' @param upper Vector of upper bounds. Upper bounds are included
#' @param probs Cell probabilities
#' @param target.prob The probability of the event, at the output sample size.
#' @param method Method used for computation. Only method currently implemented is "exact"
#' @return The sample size parameter at which the the target probability of the given event is achieved.
#'
#' @details
#' If only lower is given, then the result is the smallest size such that
#' pmultinom(lower=lower, size=size, probs=probs) >= target.prob. If only upper
#' is given, then the result is the smallest size such that
#' pmultinom(upper=upper, size=size, probs=probs) <= target.prob. Behavior when
#' both lower and upper are given is not yet implemented.
#'
#' @examples
#' # How many cells must be sequenced to have a 95% chance of
#' # observing at least 2 from each subclone of a tumor? (Data
#' # from Casasent et al (2018); see vignette("pmultinom") for
#' # details of this example)
#' 
#' # Input: 
#' ncells <- 204
#' subclone.freqs <- c(43, 20, 82, 17, 5, 37)/ncells
#' target.number <- c(2, 2, 2, 2, 2, 0)
#' lower.bound <- target.number - 1
#' invert.pmultinom(lower=lower.bound, probs=subclone.freqs,
#'                  target.prob=.95, method="exact")
#' # Output:
#' # [1] 192
#' 
#' @references 
#'     Casasent, A. K., Schalck, A., Gao, R., Sei, E., Long, A., Pangburn, W., ... & Navin, N. E. (2018). Multiclonal Invasion in Breast Tumors Identified by Topographic Single Cell Sequencing. Cell. doi:10.1016/j.cell.2017.12.007
#' 
#' @seealso \code{\link{pmultinom}}
#'
#' @export
invert.pmultinom <- function(lower=-Inf, upper=Inf, probs, target.prob, method)
{
  if (!isTRUE(method %in% c("exact", "edgeworth"))) stop('Method must be specified. Only available method right now is "exact"')
  # Checking input
  # Check that either the input is a numeric vector (which may include NA values), or it is entirely NA:
  stopifnot(is.numeric(lower) | all(is.na(lower)))
  stopifnot(is.numeric(upper) | all(is.na(upper)))
  stopifnot(is.numeric(probs) | all(is.na(probs)))
  stopifnot(is.numeric(target.prob) | all(is.na(target.prob)))
  # Check that the non-missing probabilities are really probabilities:
  stopifnot(all(sapply(0 <= probs & probs <= 1, isTRUE) | is.na(probs)))
  stopifnot(all(sapply(0 <= target.prob & target.prob <= 1, isTRUE) | is.na(target.prob)))
  # Test that the arguments recycle properly
  tryCatch(mapply(function(x, y, z) NULL, lower, upper, probs),
           warning=function(w) stop(sprintf("Error generated from warning in %s: %s",
                                            format(conditionCall(w)), conditionMessage(w))))
  # Return NA for NA input
  if (all(is.na(target.prob)) | any(is.na(c(lower, upper, probs))))
  {
    return(NA)
  }
  # The inverse value may not exist if both lower and upper are given. For the
  # time being, rather than determining existence or nonexistence, I am just
  # leaving the function undefined for this case, and raising an error when both
  # are given.
  if (!isTRUE((all(is.infinite(lower) & lower < 0) | all(is.infinite(upper) & upper > 0))))
  {
    stop("Lower and upper both given. Inverse is not guaranteed to exist under these conditions, and correct behavior has not yet been implemented.")
  }
  # If neither are given, then the answer is 0.
  # Note that this condition will fail if there are any NA's in lower and upper.
  if (isTRUE(all(is.infinite(lower) & lower < 0) & all(is.infinite(upper) & upper > 0)))
  {
    return(0)
  }
  # Exactly of the following two conditions has to be true for us to have passed the
  # last two cases, so the variable "increasing" is guaranteed to be defined
  # Only upper is given:
  if (all(is.infinite(lower) & lower < 0)) increasing <- FALSE
  # Only lower is given:
  if (all(is.infinite(upper) & upper > 0)) increasing <- TRUE
  # Make recycled versions of the arguments
  r.lower <- mapply(function(x, y, z) x, lower, upper, probs)
  r.upper <- mapply(function(x, y, z) y, lower, upper, probs)
  r.probs <- mapply(function(x, y, z) z, lower, upper, probs)
  # Check that the probabilities sum to 1, after recycling
  if (!isTRUE(all.equal(sum(r.probs), 1))) stop(sprintf("Probabilities do not sum to 1 after recycling: %s", paste(r.probs, collapse=", ")))
  # Infinite loops are possible in this function of the condition is impossible
  # to satisfy. Better catch that early and raise an error.
  # The only infinite loop I can think of is a lower bound on a zero probability category.
  is.zero.prob <- sapply(r.probs, function(p) isTRUE(all.equal(p, 0)))
  is.lower.bound <- r.lower >= 0
  causes.infinite.loop <- sapply(is.zero.prob & is.lower.bound, isTRUE)
  if (any(causes.infinite.loop))
  {
    stop(sprintf("Requiring more than %.3f in a category with probability %.3f will cause an infinite loop", r.lower[causes.infinite.loop][1], r.probs[causes.infinite.loop][1]))
  }

  if (method == "exact")
  {
    # How far to skip ahead on each iteration
    a <- 3
    # Generate values up to n=a^2
    v <- a
    n <- v^2
    distribution <- sum.tpois.pmf(r.lower, r.upper, n * r.probs, n)
    path <- prob.between(r.lower, r.upper, r.probs, 0:n, distribution, n)
    while ((increasing & path[length(path)] < max(target.prob, na.rm=TRUE)) |
           !increasing & path[length(path)] > min(target.prob, na.rm=TRUE))
    {
      # Increment n, repeat
      old.n <- n
      v <- v + a
      n <- v^2 - 1
      distribution <- sum.tpois.pmf(r.lower, r.upper, n * r.probs, n)
      # Calculate the probability of exceeding bs
      #path <- c(path, prob.greater(bs, freqs, (old.n+1):n, distribution, n))
      path <- c(path, prob.between(r.lower, r.upper, r.probs, (old.n+1):n, distribution, n))
    }
    # The index of the path where the condition is first satisfied
    index.condition.satisfied <-
      sapply(target.prob, function(tp)
        which(increasing & path >= tp |
             !increasing & path <= tp)[1])
    # R is one-indexed, so an index of 1 represents a sample size of 0
    sample.size.condition.satisfied <- index.condition.satisfied - 1
    return(sample.size.condition.satisfied)
    } else if (method == "edgeworth")
    {
      if (!all(upper == Inf)) stop("Edgeworth expansion only implemented for lower bounds so far")
      # Find a starting value by inverting an approximation If you try
      # putting negative numbers in this it will crash your R session, so
      # the input has to be filtered, which shouldn't change the answer
      # since those conditions are automatically satisfied anyway
      starting.value <-
          invert.poisson.survival.function(probs[lower>=0], lower[lower>=0], rep.int(1, length(probs))[lower>=0], target.prob, eps=10^-2)
      # Newton-raphson rootfinding
      # This uses way more function values than necessary, so either hte
      # algorithm should be modified, or the closure should be memoized
      current.in <- round(starting.value)
      prob.closure <- function(l) pmultinom(lower=lower, probs=probs, size=l, method="edgeworth")
      while (TRUE)
      {
        current.out <- prob.closure(current.in) - target.prob
        below <- prob.closure(current.in - 1) - target.prob
        above <- prob.closure(current.in + 1) - target.prob
        
        if ((current.out > 0 & below < 0)
            | current.out == 0) break
        
        derivative <- (above - below) / 2
        
        next.in <- round(current.in - current.out / derivative)
        if (next.in == current.in)
        {
          if (current.out < 0)
          {
            current.in <- current.in + 1
          } else
          {
            current.in <- current.in - 1
          }
        } else
        {
          current.in <- next.in
        }
      }
      return(current.in)
    }
}
