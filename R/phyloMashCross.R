#' @title Cross-species shift versus within-species condition effect, with mash
#'
#' @description
#' One focal species has samples in two conditions (e.g. case and control);
#' the other species are "controls" of the cross-species comparison. For every
#' gene, two effects are estimated jointly in one phylogenetic linear model
#' with precision weights (\code{\link{phyloWeightedFit}}):
#' \itemize{
#' \item \eqn{\beta_{evo}}: the shift on the terminal branch of the focal species,
#' i.e. the deviation of its control samples from the phylogenetic expectation
#' given the other species;
#' \item \eqn{\beta_{case}}: the condition effect within the focal species
#' (case minus control).
#' }
#' Then \code{\link[mashr]{mash}} on the genes by 2 matrix of effects, with the
#' exact sampling correlation of the two estimates, and the canonical patterns
#' of a two-condition analysis (as \code{dreamlet::run_mash}): null; evolutionary
#' shift only; condition effect only; both, same sign (the condition recapitulates
#' the evolutionary shift); both, opposite sign (the condition reverts it); both,
#' independent. These patterns form the mash prior. Genes are then classified
#' from the posterior of the two effects: with \eqn{p_e = 1 - lfsr_{evo}} and
#' \eqn{p_c = 1 - lfsr_{case}} the posterior probabilities that each effect is
#' non-zero with its posterior-mean sign, \code{class_prob} is
#' \eqn{(1-p_e)(1-p_c)} for \code{null}, \eqn{p_e (1 - p_c)} for \code{evo_only},
#' \eqn{(1-p_e) p_c} for \code{case_only}, \eqn{p_e p_c} for \code{shared_same}
#' or \code{shared_opposite} according to the posterior-mean signs, where
#' \eqn{p^0} is the posterior probability that the effect is exactly zero (weight of
#' the components with zero variance for that coordinate); the probabilities are
#' normalised to sum to one (the sign posteriors are treated as independent given
#' the fit; \code{independent}, if used, is only a prior pattern and gets no class probability). The pattern-level
#' posterior weights are also returned (\code{pattern_prob}).
#'
#' @param counts genes by samples matrix of counts.
#' @param phy the species tree (one tip per species), ultrametric.
#' @param col_species the species of each sample.
#' @param case a 0/1 vector (one per sample): 1 for the case samples of the focal
#' species, 0 for all other samples (controls of the focal species and all the
#' samples of the other species).
#' @param focal the focal species (the only species with \code{case == 1}).
#' @param model \code{"BM"} (the default) or \code{"OUfixedRoot"}.
#' @param other_branches if \code{TRUE}, the shifts on the other branches are also
#' estimated (as in \code{\link{phyloWeightedFit}}) and returned in \code{fit}.
#' @param shared_ratios variance ratios var(case) / var(evo) of the "both" patterns.
#' Default \code{c(0.5, 1, 2)}.
#' @param shared_cor correlations of the "both" patterns; \code{1} (the default,
#' rank-one patterns: the two effects are proportional). Adding intermediate values
#' (e.g. \code{0.5}) or an "independent" pattern (\code{independent = TRUE}) makes
#' the shared classes absorb genes with a single effect and noise in the other,
#' so that the \code{evo_only} and \code{case_only} classes become uncallable;
#' see the package benchmarks.
#' @param independent add a pattern with both effects non-zero and uncorrelated. Default \code{FALSE}.
#' @param n_grid number of mash scales.
#' @param ... further arguments passed to \code{\link{phyloWeightedFit}}.
#'
#' @return An object of class \code{PhyloMashCross}, a list with the \code{mash}
#' fit, the genes by 2 matrices \code{Bhat}, \code{Shat}, \code{lfsr}, \code{pm},
#' the genes by classes matrix \code{class_prob} (columns \code{null},
#' \code{evo_only}, \code{case_only}, \code{shared_same}, \code{shared_opposite},
#' \code{independent}), \code{top_class}, \code{pattern_weights}, the sampling
#' correlation \code{V}, and the underlying weighted fit \code{fit}.
#'
#' @seealso \code{\link{phyloWeightedFit}}, \code{\link{phyloMash}}, \code{\link{classifyCross}}
#'
#' @export
#'
phyloMashCross <- function(counts, phy, col_species, case, focal, model = c("BM", "OUfixedRoot"),
                           other_branches = FALSE, shared_ratios = c(0.5, 1, 2), shared_cor = 1, independent = FALSE,
                           n_grid = 15, ...) {
  if (!requireNamespace("mashr", quietly = TRUE)) stop("Package 'mashr' is needed for function 'phyloMashCross'.")
  model <- match.arg(model)
  case <- as.numeric(case)
  if (any(case == 1 & col_species != focal)) stop("`case == 1` is only allowed for samples of the focal species.")
  if (sum(case == 1) < 2 || sum(case == 0 & col_species == focal) < 2) stop("The focal species needs at least two case and two control samples.")
  n <- ncol(counts)
  ## design: intercept + case; the focal terminal branch indicator is the tested "shift" (controls + cases both descend from it,
  ## so the branch shift is identified from the controls given the case covariate)
  design <- cbind(1, case); dimnames(design) <- list(colnames(counts), c("(Intercept)", "case"))
  fit <- phyloWeightedFit(counts, phy, col_species, design = design, shifts = TRUE, model = model, ...)
  if (!focal %in% colnames(fit$Bhat)) stop("The focal species' terminal branch is not testable (check the tree).")
  ## joint estimate of (evo, case): refit the pair at the posterior-mean lambda of the focal-branch test, for the exact
  ## sampling correlation. Both coefficients in one design: intercept + focal branch + case.
  cor_pair <- cross_sampling_correlation(fit, focal)
  Bhat <- cbind(evo = fit$Bhat[, focal], case = fit$Bhat[, "case"])
  Shat <- cbind(evo = fit$Shat[, focal], case = fit$Shat[, "case"])
  V <- matrix(c(1, cor_pair, cor_pair, 1), 2, dimnames = list(c("evo", "case"), c("evo", "case")))
  data <- mashr::mash_set_data(Bhat, Shat, df = fit$df, V = V)
  ## patterns: singletons; and for each variance ratio r = var(case)/var(evo), correlations
  ## +1 / +0.5 (same sign), -1 / -0.5 (opposite sign) and 0 (independent)
  Ulist <- list(evo_only = matrix(c(1, 0, 0, 0), 2), case_only = matrix(c(0, 0, 0, 1), 2))
  tag <- c(evo_only = "evo_only", case_only = "case_only")
  for (r in shared_ratios) {
    for (rho in shared_cor) {
      nm <- paste0("shared_same_r", r, "_c", rho); Ulist[[nm]] <- matrix(c(1, rho * sqrt(r), rho * sqrt(r), r), 2); tag[nm] <- "shared_same"
      nm <- paste0("shared_opposite_r", r, "_c", rho); Ulist[[nm]] <- matrix(c(1, -rho * sqrt(r), -rho * sqrt(r), r), 2); tag[nm] <- "shared_opposite"
    }
    if (independent) { Ulist[[paste0("independent_r", r)]] <- matrix(c(1, 0, 0, r), 2); tag[paste0("independent_r", r)] <- "independent" }
  }
  g_auto <- mashr:::autoselect_grid(data, sqrt(2))
  grid <- exp(seq(log(min(g_auto)), log(max(g_auto)), length.out = n_grid))
  m <- mashr::mash(data, Ulist, grid = grid, outputlevel = 2, verbose = FALSE)
  w <- m$posterior_weights
  base <- sub("\\.[0-9]+$", "", colnames(w))
  cls <- c("null", "evo_only", "case_only", "shared_same", "shared_opposite", "independent")
  tag_w <- ifelse(base == "null", "null", tag[base])
  pattern_prob <- sapply(cls, function(k) rowSums(w[, which(tag_w == k), drop = FALSE]))
  if (is.null(dim(pattern_prob))) pattern_prob <- matrix(pattern_prob, nrow = nrow(Bhat))
  dimnames(pattern_prob) <- list(rownames(Bhat), cls)
  pi <- m$fitted_g$pi
  pattern_weights <- tapply(pi, ifelse(names(pi) == "null", "null", tag[sub("\\.[0-9]+$", "", names(pi))]), sum)
  ## class probabilities from the posterior of the effects themselves: P(evo != 0 with sign s_e, case != 0 with sign s_c)
  ## from the per-coordinate posterior sign probabilities (lfsr) and the posterior correlation of the pair
  pm <- ashr::get_pm(m); psd <- ashr::get_psd(m); lfsr <- ashr::get_lfsr(m)
  ## P(effect exactly zero) per coordinate = posterior weight on the components with zero variance for that coordinate
  ## (null, and the singleton of the other coordinate); P(nonzero with the posterior-mean sign) = 1 - lfsr - that
  zero_evo <- base %in% c("null", "case_only"); zero_case <- base %in% c("null", "evo_only")
  p0_evo <- rowSums(w[, zero_evo, drop = FALSE]); p0_case <- rowSums(w[, zero_case, drop = FALSE])
  p_evo <- 1 - lfsr[, "evo"]; p_case <- 1 - lfsr[, "case"]
  same <- sign(pm[, "evo"]) == sign(pm[, "case"])
  class_prob <- cbind(null = p0_evo * p0_case,
                      evo_only = p_evo * p0_case,
                      case_only = p0_evo * p_case,
                      shared_same = p_evo * p_case * same,
                      shared_opposite = p_evo * p_case * !same,
                      independent = 0)
  class_prob <- class_prob / pmax(rowSums(class_prob), 1e-12)
  dimnames(class_prob) <- list(rownames(Bhat), cls)
  res <- list(mash = m, Bhat = Bhat, Shat = Shat, lfsr = lfsr, pm = pm, psd = psd,
              class_prob = class_prob, pattern_prob = pattern_prob,
              top_class = cls[max.col(class_prob, ties.method = "first")],
              pattern_weights = pattern_weights, V = V, fit = fit, focal = focal, genes = rownames(Bhat))
  class(res) <- "PhyloMashCross"
  res
}

#' @title Sampling correlation of the focal-branch and case estimates
#'
#' @description
#' Under the null, the estimate of the shift on the focal branch (from the
#' design \code{intercept + case}) and the estimate of the case effect (same
#' design) are correlated because the case samples enter both. The correlation
#' is computed from the whitened design at the posterior-mean correlation
#' structure of the focal-branch test (one value for all genes, with the
#' average weight profile).
#'
#' @param fit a \code{PhyloShiftFit} from \code{\link{phyloWeightedFit}} with a \code{case} coefficient.
#' @param focal the focal species.
#'
#' @return The correlation.
#'
#' @keywords internal
#'
cross_sampling_correlation <- function(fit, focal) {
  n <- nrow(fit$Z)
  lam <- stats::median(fit$lambda[, focal])
  tr_u <- rescale_tree(fit$phy_rep)
  alpha <- if (fit$model == "OUfixedRoot") rho_prime_inv(1 - stats::median(fit$rho[, focal]), 1, getBoundsSelectionStrength(tr_u)) else 0
  C <- ape::vcv(eb_transform_tree(tr_u, fit$model, 1, alpha))
  Dn <- colMeans((1 / fit$weights) / rowMeans(1 / fit$weights))
  Sig <- lam * C + (1 - lam) * diag(Dn, n)
  R <- chol(Sig)
  X <- cbind(fit$design[, 1], fit$Z[, focal], fit$design[, "case"])
  Xt <- backsolve(R, X, transpose = TRUE)
  XtX_inv <- solve(crossprod(Xt))
  stats::cov2cor(XtX_inv)[2, 3]
}

#' @title Classify genes from a cross analysis at a Bayesian FDR
#'
#' @description
#' Call, for each class of \code{\link{phyloMashCross}} other than \code{null},
#' the genes with the highest posterior probability of that class while keeping
#' the expected proportion of wrong calls below \code{level}.
#'
#' @param object a \code{PhyloMashCross} object.
#' @param level the Bayesian false discovery rate. Default 0.05.
#'
#' @return A genes by classes logical matrix.
#'
#' @export
#'
classifyCross <- function(object, level = 0.05) {
  cp <- object$class_prob[, colnames(object$class_prob) != "null", drop = FALSE]
  out <- cp > 2
  for (k in colnames(cp)) {
    p <- cp[, k]; o <- order(-p); cum <- cumsum(1 - p[o]) / seq_along(o)
    kk <- suppressWarnings(max(which(cum <= level)))
    sel <- rep(FALSE, length(p)); if (is.finite(kk)) sel[o[seq_len(kk)]] <- TRUE
    out[, k] <- sel
  }
  out
}
