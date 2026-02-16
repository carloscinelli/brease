#' Print Method for BREASE Objects
#'
#' @description Prints a concise summary of the BREASE analysis results.
#'
#' @param x An object of class \code{"brease"}.
#' @param digits Number of significant digits to display (default 3).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.brease <- function(x, digits = 3, ...) {
  method <- x$info$method
  method_label <- switch(method,
    brease = if (x$info$mono) "BREASE (monotone)" else "BREASE",
    ib     = "Independent Beta",
    lt     = "Logit-Transform",
    method
  )

  cat("Bayesian Analysis of Binary Experiment\n")
  cat("Method:", method_label, "\n\n")

  # Data
  cat("Data:\n")
  cat(sprintf("  Control:   y0 = %d, N0 = %d (rate = %.4f)\n",
              x$data$y0, x$data$N0, x$data$y0 / x$data$N0))
  cat(sprintf("  Treatment: y1 = %d, N1 = %d (rate = %.4f)\n",
              x$data$y1, x$data$N1, x$data$y1 / x$data$N1))
  cat("\n")

  # Prior
  cat("Prior:\n")
  if (method == "brease") {
    p <- x$prior
    cat(sprintf("  theta0 ~ Beta(%.2f, %.2f)  [mu0 = %.2f, n0 = %.1f]\n",
                p$mu0 * p$n0, (1 - p$mu0) * p$n0, p$mu0, p$n0))
    cat(sprintf("  eta_e  ~ Beta(%.2f, %.2f)  [mue = %.2f, ne = %.1f]\n",
                p$mue * p$ne, (1 - p$mue) * p$ne, p$mue, p$ne))
    if (!x$info$mono) {
      cat(sprintf("  eta_s  ~ Beta(%.2f, %.2f)  [mus = %.2f, ns = %.1f]\n",
                  p$mus * p$ns, (1 - p$mus) * p$ns, p$mus, p$ns))
    } else {
      cat("  eta_s  = 0 (monotone/no-harm)\n")
    }
  } else if (method == "ib") {
    p <- x$prior
    cat(sprintf("  theta0 ~ Beta(%.2f, %.2f)\n", p$a[1], p$b[1]))
    cat(sprintf("  theta1 ~ Beta(%.2f, %.2f)\n", p$a[2], p$b[2]))
  } else if (method == "lt") {
    p <- x$prior
    cat(sprintf("  psi  ~ N(%.1f, %.1f)\n", p$mu_psi, p$sigma_psi))
    cat(sprintf("  beta ~ N(%.1f, %.1f)\n", p$mu_beta, p$sigma_beta))
  }
  cat("\n")

  # Posterior summaries (key quantities)
  cat("Posterior Summaries:\n")
  ps <- x$post.summaries
  display_rows <- c("theta0", "theta1", "rr", "ef", "or", "rd")
  display_rows <- display_rows[display_rows %in% rownames(ps)]
  labels <- c(theta0 = "  Baseline Risk (theta0)",
              theta1 = "  Treatment Risk (theta1)",
              rr     = "  Risk Ratio (RR)",
              ef     = "  Efficacy Fraction (EF)",
              or     = "  Odds Ratio (OR)",
              rd     = "  Risk Difference (RD)")
  for (row in display_rows) {
    vals <- ps[row, ]
    cat(sprintf("%-30s  Mean: %s  [%s, %s]\n",
                labels[row],
                formatC(vals["mean"], format = "g", digits = digits),
                formatC(vals["lw"], format = "g", digits = digits),
                formatC(vals["up"], format = "g", digits = digits)))
  }
  cat("\n")

  # Bayes Factor
  cat(sprintf("Bayes Factor (BF10): %s\n",
              formatC(x$bayes.factors$bf10, format = "g", digits = digits)))
  cat(sprintf("Bayes Factor (BF01): %s\n",
              formatC(x$bayes.factors$bf01, format = "g", digits = digits)))

  invisible(x)
}


#' Summary Method for BREASE Objects
#'
#' @description Produces a detailed summary of the BREASE analysis results
#'   including all posterior quantities.
#'
#' @param object An object of class \code{"brease"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the full posterior summaries matrix.
#'
#' @export
summary.brease <- function(object, ...) {
  cat("Posterior Summaries:\n\n")
  ps <- object$post.summaries
  print(round(ps, 5))
  cat("\n")
  cat(sprintf("Bayes Factor BF10 = %.4f\n", object$bayes.factors$bf10))
  cat(sprintf("Bayes Factor BF01 = %.4f\n", object$bayes.factors$bf01))
  invisible(ps)
}
