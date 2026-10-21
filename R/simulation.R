#' Simulate data from a forward continuation-ratio model
#'
#' Generates a predictor matrix and an ordinal outcome under the forward
#' continuation-ratio (FCR) model used in the fusedFCR simulation study.
#' Two designs are:
#' \itemize{
#'   \item \code{time_varying = FALSE} (non-temporal): predictors are i.i.d.
#'   standard normal (column-standardized), and the active genes have constant
#'   or step-shaped effects across transitions.
#'   \item \code{time_varying = TRUE} (temporal): cells are ordered along a
#'   latent time in [0, 1]; genes 1-3 are early, transient, and late markers
#'   whose expression follows smooth trajectories, and their effects change
#'   sign across transitions.
#' }
#' The outcome is drawn sequentially: at transition \eqn{t}, a cell still at
#' risk moves to stage \eqn{t} with probability
#' \eqn{\mathrm{logit}^{-1}(\alpha_t + x^\top \beta_t)}; cells that never
#' transition end in stage \code{T_max}.
#'
#' @param N Integer number of cells (observations).
#' @param P Integer number of genes (predictors). Must be at least \code{n_nonz}.
#' @param T_max Integer number of ordinal stages. Only \code{T_max = 6} (five
#'   transitions) is currently implemented.
#' @param n_nonz Integer number of truly active genes. For
#'   \code{time_varying = FALSE}, one of 3, 6, or 9; for
#'   \code{time_varying = TRUE}, only 3.
#' @param time_varying Logical; if \code{TRUE}, simulate the temporal design,
#'   otherwise the non-temporal design. Default \code{FALSE}.
#' @param seed Integer random seed. Default 2026.
#'
#' @return A list with elements
#'   \item{X}{N x P predictor matrix with columns \code{Gene_1}, ..., \code{Gene_P}.}
#'   \item{Y}{Integer vector of length N with stages in \code{1:T_max}.}
#'   \item{true_beta}{P x (T_max - 1) matrix of true transition-specific coefficients.}
#'   \item{alpha}{Length (T_max - 1) vector of true intercepts.}
#'   \item{time_varying, n_nonz, seed}{The settings used.}
#'
#' @examples
#' sim <- generate_sim_data(N = 500, P = 20, time_varying = FALSE)
#' table(sim$Y)
#'
#' sim_t <- generate_sim_data(N = 500, P = 20, time_varying = TRUE)
#' table(sim_t$Y)
#'
#' @export
generate_sim_data <- function(N, P, T_max = 6, n_nonz = 3, time_varying = FALSE, seed = 2026) {
  if (!is.logical(time_varying) || length(time_varying) != 1 || is.na(time_varying)) {
    stop("time_varying must be TRUE or FALSE.")
  }
  if (T_max != 6) {
    stop("Only T_max = 6 (five transitions) is currently implemented.")
  }
  if (P < n_nonz) {
    stop("P must be at least n_nonz.")
  }
  TT <- T_max - 1
  set.seed(seed)
  if (time_varying) {
    sim <- sim_temporal(N = N, P = P, TT = TT, n_nonz = n_nonz)
  } else {
    sim <- sim_nontemporal(N = N, P = P, TT = TT, n_nonz = n_nonz)
  }
  X <- sim$X
  true_beta <- sim$true_beta
  alpha <- sim$alpha
  # Sequential FCR draw: a cell still at risk moves to stage t with prob logit^-1(alpha_t + x'beta_t)
  Y <- rep(T_max, N)
  alive <- rep(TRUE, N)
  for (t in 1:TT) {
    eta_t <- alpha[t] + X %*% true_beta[, t]
    prob_t <- 1 / (1 + exp(-eta_t))
    transitioned_now <- stats::rbinom(N, size = 1, prob = prob_t) == 1
    event_occured <- alive & transitioned_now
    Y[event_occured] <- t
    alive[event_occured] <- FALSE
  }
  list(X = X, Y = Y, true_beta = true_beta, alpha = alpha, time_varying = time_varying, n_nonz = n_nonz, seed = seed)
}

#' Non-temporal simulation design
#'
#' i.i.d. standard normal predictors, column-standardized, with constant or
#' step-shaped effects for the active genes.
#'
#' @param N Number of cells.
#' @param P Number of genes.
#' @param TT Number of transitions (T_max - 1).
#' @param n_nonz Number of active genes: 3, 6, or 9.
#' @return A list with X, true_beta, and alpha.
#' @keywords internal
#' @noRd
sim_nontemporal <- function(N, P, TT, n_nonz) {
  X <- matrix(stats::rnorm(N * P), nrow = N, ncol = P)
  X <- scale(X)
  colnames(X) <- paste0("Gene_", 1:P)
  true_beta <- matrix(0, nrow = P, ncol = TT)
  if (n_nonz == 3) {
    true_beta[1, ] <- c(1.0, 1.0, 1.0, 1.0, 1.0)
    true_beta[2, ] <- c(0.0, 0.0, 1.0, 1.0, 1.0)
    true_beta[3, ] <- c(-1.0, -1.0, -1.0, -1.0, -1.0)
    alpha <- c(-2.0, -1.5, -1.0, -0.5, 1.2)
  } else if (n_nonz == 6) {
    true_beta[1, ] <- c(1.0, 1.0, 1.0, 1.0, 1.0)
    true_beta[2, ] <- c(0.0, 0.0, 1.0, 1.0, 1.0)
    true_beta[3, ] <- c(-1.0, -1.0, -1.0, -1.0, -1.0)
    for (bb in 4:6) {
      true_beta[bb, ] <- c(0.0, 0.0, 0.0, 1.0, 1.0)
    }
    alpha <- c(-2.0, -1.5, -1.2, 0.0, 2.0)
  } else if (n_nonz == 9) {
    true_beta[1, ] <- c(1.0, 1.0, 1.0, 1.0, 1.0)
    true_beta[2, ] <- c(0.0, 0.0, 2.0, 2.0, 2.0)
    true_beta[3, ] <- c(-1.0, -1.0, -1.0, -1.0, -1.0)
    for (bb in 4:8) {
      true_beta[bb, ] <- c(0.0, 0.0, 0.0, 1.0, 1.0)
    }
    true_beta[9, ] <- c(-1.0, -1.0, -1.0, -1.0, -1.0)
    alpha <- c(-2.3, -1.5, -1.5, 0.0, 3.0)
  } else {
    stop("For time_varying = FALSE, n_nonz must be 3, 6, or 9.")
  }
  list(X = X, true_beta = true_beta, alpha = alpha)
}

#' Temporal simulation design
#'
#' Cells are ordered along a latent time in [0, 1]. Gene 1 is an early marker
#' (exponential decay), gene 2 a transient marker (Gaussian bump at mid-time),
#' and gene 3 a late marker (exponential rise); all other genes are noise.
#' Effects change sign across transitions so each marker drives the
#' transitions at its own stage.
#'
#' @param N Number of cells.
#' @param P Number of genes.
#' @param TT Number of transitions (T_max - 1).
#' @param n_nonz Number of active genes: only 3.
#' @return A list with X, true_beta, and alpha.
#' @keywords internal
#' @noRd
sim_temporal <- function(N, P, TT, n_nonz) {
  if (n_nonz != 3) {
    stop("For time_varying = TRUE, only n_nonz = 3 is implemented.")
  }
  latent_time <- seq(0, 1, length.out = N)
  X <- matrix(stats::rnorm(N * P, mean = 0, sd = 0.20), nrow = N, ncol = P)
  # Early marker
  X[, 1] <- X[, 1] + 4.5 * exp(-5.5 * latent_time)
  # Transient (mid) marker
  X[, 2] <- X[, 2] + 5.0 * exp(-((latent_time - 0.5)^2) / (2 * 0.11^2))
  # Late marker
  X[, 3] <- X[, 3] + 4.5 * exp(-5.5 * (1 - latent_time))
  colnames(X) <- paste0("Gene_", 1:P)
  true_beta <- matrix(0, nrow = P, ncol = TT)
  true_beta[1, ] <- c(5.5, 3.0, 0.0, -3.0, -5.5)
  true_beta[2, ] <- c(-2.5, 2.0, 5.5, 2.0, -2.5)
  true_beta[3, ] <- c(-5.5, -3.0, 0.0, 3.0, 5.5)
  alpha <- c(-4.5, -4.2, -4.0, -3.2, -2.0)
  list(X = X, true_beta = true_beta, alpha = alpha)
}