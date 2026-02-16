#' Aspirin Myocardial Infarction Trial
#'
#' @description Data from the Physicians' Health Study (1989), a large-scale
#'   randomized trial studying the effect of aspirin on myocardial infarction
#'   (MI). 22,071 male physicians were randomized to aspirin or placebo.
#'
#' @format A data frame with 1 row and 4 variables:
#' \describe{
#'   \item{y0}{Number of MIs in the placebo group (26).}
#'   \item{N0}{Number of subjects in the placebo group (11,034).}
#'   \item{y1}{Number of MIs in the aspirin group (10).}
#'   \item{N1}{Number of subjects in the aspirin group (11,037).}
#' }
#'
#' @references
#' Steering Committee of the Physicians' Health Study Research Group (1989).
#' Final report on the aspirin component of the ongoing Physicians' Health
#' Study. \emph{New England Journal of Medicine}, 321(3), 129--135.
#'
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' data(aspirin)
#' brease(aspirin$y0, aspirin$y1, aspirin$N0, aspirin$N1)
"aspirin"


#' NEJM Clinical Trials Dataset
#'
#' @description A collection of 39 clinical trials published in the
#'   \emph{New England Journal of Medicine}, used to compare Bayesian
#'   analysis methods. Each row is a separate trial with binary outcomes.
#'
#' @format A data frame with 39 rows and 10 variables:
#' \describe{
#'   \item{No.}{Trial reference number.}
#'   \item{P.value}{Frequentist p-value.}
#'   \item{BF01}{Bayes factor in favor of H0.}
#'   \item{y1}{Number of events in group 1 (control).}
#'   \item{y2}{Number of events in group 2 (treatment).}
#'   \item{n1}{Number of subjects in group 1 (control).}
#'   \item{n2}{Number of subjects in group 2 (treatment).}
#'   \item{y1n1}{Event rate in group 1 (y1/n1).}
#'   \item{y2n2}{Event rate in group 2 (y2/n2).}
#'   \item{ynavg}{Average event rate across both groups.}
#' }
#'
#' @references
#' Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
#' Experiments. \emph{Bayesian Analysis}. \doi{10.1214/25-BA1506}
#'
#' @examples
#' data(nejm)
#' head(nejm)
"nejm"
