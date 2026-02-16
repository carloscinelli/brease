#' Bounds on Probabilities of Causation
#'
#' @description Computes nonparametric bounds on the probabilities of
#'   causation (efficacy and side-effect risk) from the observed data,
#'   based on plug-in estimates of \eqn{\theta_0} and \eqn{\theta_1}.
#'
#' @param y0 Number of events in the control group.
#' @param y1 Number of events in the treatment group.
#' @param N0 Number of subjects in the control group.
#' @param N1 Number of subjects in the treatment group.
#'
#' @return A list with components:
#' \describe{
#'   \item{eta_s_lower}{Lower bound on the risk of side effects.}
#'   \item{eta_s_upper}{Upper bound on the risk of side effects.}
#'   \item{eta_s_mid}{Midpoint of the side-effect risk bounds.}
#'   \item{eta_e_lower}{Lower bound on efficacy.}
#'   \item{eta_e_upper}{Upper bound on efficacy.}
#'   \item{eta_e_mid}{Midpoint of the efficacy bounds.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' # Aspirin trial
#' brease_bounds(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
#'
#' # COVID vaccine trial
#' brease_bounds(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
#'
#' @export
brease_bounds <- function(y0, y1, N0, N1) {
  check_data(y0, y1, N0, N1)

  # Add 1 success and 1 failure for numerical stability
  theta0_hat <- (y0 + 1) / (N0 + 2)
  theta1_hat <- (y1 + 1) / (N1 + 2)

  # Bounds on eta_s (side-effect risk)
  if (theta0_hat == 1) {
    eta_s_lower <- 0
    eta_s_upper <- 1
  } else {
    eta_s_lower <- max(0, (theta1_hat - theta0_hat) / (1 - theta0_hat))
    eta_s_upper <- min(1, theta1_hat / (1 - theta0_hat))
  }
  eta_s_mid <- (eta_s_lower + eta_s_upper) / 2

  # Bounds on eta_e (efficacy)
  if (theta0_hat == 0) {
    eta_e_lower <- 0
    eta_e_upper <- 1
  } else {
    eta_e_lower <- max(0, (theta0_hat - theta1_hat) / theta0_hat)
    eta_e_upper <- min(1, (1 - theta1_hat) / theta0_hat)
  }
  eta_e_mid <- (eta_e_lower + eta_e_upper) / 2

  list(
    eta_s_lower = eta_s_lower,
    eta_s_upper = eta_s_upper,
    eta_s_mid   = eta_s_mid,
    eta_e_lower = eta_e_lower,
    eta_e_upper = eta_e_upper,
    eta_e_mid   = eta_e_mid
  )
}
