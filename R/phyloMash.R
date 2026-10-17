#' @title Phylogenetic shift detection on species means with mash
#'
#' @description
#' Detect branch-specific shifts in expression (or accessibility) across
#' species with replicates, using \code{\link[mashr]{mash}} on the species means.
#' The species means of a gene are the sufficient statistics of its species-level
#' effects; their deviations from the across-species mean are the "effects" of a
#' mash analysis, and the prior on these effects is a mixture of scaled
#' covariance patterns learned from all genes:
#' \itemize{
#' \item "drift" patterns: no shift, the deviations follow the phylogenetic
#' process (BM, or OU for several values of \eqn{\rho' = 1 - \rho}, see
#' \code{\link{rhoFromAlpha}}), possibly plus species-specific
#' non-phylogenetic variance (\code{noise_ratios}) or only such variance
#' (\code{identity});
#' \item "shift" patterns: BM drift plus a shift shared by all the species
#' descending from one branch, with a shift variance equal to
#' \code{shift_ratios} times the drift variance, or a pure shift.
#' }
#' The within-species variance (replicates) gives the standard errors of the
#' species means, moderated across genes with \code{\link[limma]{squeezeVar}}.
#' The scale of the drift (the gene-wise phylogenetic variance) is integrated
#' over the mash scale grid, so that no gene-wise phylogenetic parameter has to
#' be estimated: the posterior weights over the patterns and scales are the
#' empirical Bayes posterior on the gene-wise phylogenetic parameters.
#'
#' A gene carries a shift on branch \eqn{b} with posterior probability equal to
#' the posterior weight of the shift patterns of \eqn{b} with a shift-to-drift
#' ratio of at least \code{min_shift_ratio} (\code{branch_prob}); patterns with a
#' smaller ratio are indistinguishable from drift and count as "no shift".
#' See \code{\link{callShifts}} to call shifts at a given Bayesian false
#' discovery rate.
#'
#' @inheritParams phylolmFit
#' @param phy the species tree (one tip per species), dated and ultrametric.
#' @param model the phylogenetic model of the drift: \code{"BM"} (the default) or
#' \code{"OUfixedRoot"}.
#' @param rho_prime_grid for the OU model, the values of \eqn{\rho'} of the OU
#' drift patterns (the BM drift pattern is always included).
#' @param shift_ratios ratios of the shift variance to the drift variance (at the tips)
#' of the shift patterns.
#' @param min_shift_ratio shift patterns with a ratio below this value count as "no shift".
#' @param noise_ratios ratios of the species-specific non-phylogenetic variance
#' to the drift variance of the drift-plus-noise patterns.
#' @param identity if \code{TRUE} (the default), add a pure species-specific variance pattern.
#' @param trend if \code{TRUE} (the default), the moderation of the within-species
#' variances allows for a mean-variance trend (\code{\link[limma]{squeezeVar}} with \code{covariate}).
#' @param n_grid number of scales in the mash grid (log-spaced over the range chosen by \code{mash}).
#' @param fit_subset if not \code{NULL}, the number of (randomly chosen) genes used to
#' estimate the mixture weights, which are then applied to all genes
#' (faster for large data sets).
#' @param ... further arguments passed to \code{\link[mashr]{mash}}.
#'
#' @return An object of class \code{PhyloMash}, a list with
#' \describe{
#' \item{mash}{the \code{mash} fit.}
#' \item{branch_prob}{genes by testable branches matrix of posterior probabilities of a shift.}
#' \item{shift_prob}{posterior probability of a shift on any branch.}
#' \item{top_branch}{the branch with the highest posterior probability.}
#' \item{pm_species}{genes by species matrix of posterior mean deviations from the across-species mean.}
#' \item{species_means, se_species}{the species means and their (moderated) standard errors.}
#' \item{pattern_weights}{the estimated prior weight of each pattern (summed over scales).}
#' \item{branches}{the branch table (testable branches).}
#' }
#'
#' @seealso \code{\link{callShifts}}, \code{\link{branchIndicators}}, \code{\link[mashr]{mash}}
#'
#' @export
#'
phyloMash <- function(object, phy, col_species = NULL, model = c("BM", "OUfixedRoot"),
                      rho_prime_grid = c(0.1, 0.5), shift_ratios = c(1, 4, 16, 64), min_shift_ratio = 4,
                      noise_ratios = c(0.25, 1, 4), identity = TRUE, trend = TRUE, n_grid = 15,
                      fit_subset = NULL, ...) {
  if (!requireNamespace("mashr", quietly = TRUE)) stop("Package 'mashr' is needed for function 'phyloMash'.")
  model <- match.arg(model)
  y <- check_expression_matrix(object)
  Y <- y$exprs
  if (is.null(col_species)) col_species <- parse_species(phy, colnames(Y))
  if (!all(col_species %in% phy$tip.label)) stop("Some species in `col_species` are not in the tree.")
  if (!ape::is.ultrametric(phy)) stop("The tree must be ultrametric.")
  sp <- factor(col_species, levels = phy$tip.label)
  if (any(table(sp) == 0)) stop("All the species of the tree must have samples.")
  S <- nlevels(sp); r <- as.vector(table(sp)); n <- ncol(Y); G <- nrow(Y)
  if (n - S < 2) stop("At least two species must have replicates to estimate the within-species variance.")
  if (any(r == 1)) message("Species with a single sample: ", paste(levels(sp)[r == 1], collapse = ", "), ". Their within-species variance is taken from the other species (pooled).")

  ## sufficient statistics
  M <- t(apply(Y, 1, function(v) tapply(v, sp, mean)))
  colnames(M) <- levels(sp)
  resid <- Y - M[, as.integer(sp), drop = FALSE]
  s2w <- rowSums(resid^2) / (n - S)
  sq <- if (trend) limma::squeezeVar(s2w, df = n - S, covariate = rowMeans(Y)) else limma::squeezeVar(s2w, df = n - S)
  s2post <- sq$var.post
  df <- n - S + if (is.finite(sq$df.prior)) sq$df.prior else 0

  ## contrasts: deviations from the unweighted across-species mean
  L <- mashr::contrast_matrix(S, "mean", levels(sp))
  CL <- L %*% diag(1 / r) %*% t(L)
  Bhat <- M %*% t(L)
  Shat <- sqrt(s2post) %*% t(sqrt(diag(CL)))
  V <- stats::cov2cor(CL)
  data <- mashr::mash_set_data(Bhat, Shat, df = df, V = V)

  ## covariance patterns in species space, mapped to the contrast space
  bi <- branchIndicators(phy, levels(sp))
  keep <- bi$branches$testable
  Zs <- bi$Z[, keep, drop = FALSE]
  branches <- bi$branches[keep, ]
  labels <- branches$label
  tr_u <- rescale_tree(phy)
  drifts <- list(BM = ape::vcv(tr_u))
  if (model == "OUfixedRoot") {
    ab <- getBoundsSelectionStrength(tr_u)
    for (rp in rho_prime_grid) {
      al <- rho_prime_inv(rp, 1, ab)
      drifts[[paste0("OU", rp)]] <- ape::vcv(rescale_tree(phylolm::transf.branch.lengths(tr_u, "OUfixedRoot", parameters = list(alpha = al))$tree))
    }
  }
  toc <- function(U) { Uc <- L %*% U %*% t(L); Uc <- (Uc + t(Uc)) / 2; Uc / max(diag(Uc)) }
  Ulist <- list(); tag <- character()
  for (dn in names(drifts)) {
    Ulist[[paste0("drift_", dn)]] <- toc(drifts[[dn]]); tag[paste0("drift_", dn)] <- "none"
    for (q in noise_ratios) {
      nm <- paste0("drift_", dn, "_noise", q)
      Ulist[[nm]] <- toc(drifts[[dn]] + q * diag(S)); tag[nm] <- "none"
    }
  }
  if (identity) { Ulist[["identity"]] <- toc(diag(S)); tag["identity"] <- "none" }
  for (b in seq_along(labels)) {
    Ub <- tcrossprod(Zs[, b])
    Ulist[[paste0("shift_", labels[b])]] <- toc(Ub); tag[paste0("shift_", labels[b])] <- labels[b]
    for (rr in shift_ratios) {
      nm <- paste0("shift_", labels[b], "_r", rr)
      Ulist[[nm]] <- toc(drifts[["BM"]] + rr * Ub)
      tag[nm] <- if (rr >= min_shift_ratio) labels[b] else "none"
    }
  }

  ## mash fit
  g_auto <- mashr:::autoselect_grid(data, sqrt(2))
  grid <- exp(seq(log(min(g_auto)), log(max(g_auto)), length.out = n_grid))
  if (!is.null(fit_subset) && fit_subset < G) {
    sub <- sort(sample.int(G, fit_subset))
    data_sub <- mashr::mash_set_data(Bhat[sub, , drop = FALSE], Shat[sub, , drop = FALSE], df = df, V = V)
    m_sub <- mashr::mash(data_sub, Ulist, grid = grid, outputlevel = 1, verbose = FALSE, ...)
    m <- mashr::mash(data, g = ashr::get_fitted_g(m_sub), fixg = TRUE, outputlevel = 2, verbose = FALSE)
  } else {
    m <- mashr::mash(data, Ulist, grid = grid, outputlevel = 2, verbose = FALSE, ...)
  }

  ## posterior probabilities of the shift patterns
  w <- m$posterior_weights
  base <- sub("\\.[0-9]+$", "", colnames(w))
  tag_w <- tag[base]; tag_w[base == "null"] <- "none"
  branch_prob <- sapply(labels, function(b) rowSums(w[, which(tag_w == b), drop = FALSE]))
  if (is.null(dim(branch_prob))) branch_prob <- matrix(branch_prob, nrow = G)
  dimnames(branch_prob) <- list(rownames(Y), labels)
  pi <- m$fitted_g$pi
  pattern_weights <- tapply(pi, sub("\\.[0-9]+$", "", names(pi)), sum)
  pmc <- ashr::get_pm(m)
  pm_species <- cbind(pmc, -rowSums(pmc)); colnames(pm_species) <- levels(sp)

  res <- list(mash = m, branch_prob = branch_prob, shift_prob = rowSums(branch_prob),
              top_branch = labels[max.col(branch_prob, ties.method = "first")],
              pm_species = pm_species, species_means = M, se_species = sqrt(s2post) %*% t(1 / sqrt(r)),
              pattern_weights = pattern_weights, branches = branches, L = L, df = df,
              model = model, phy = phy, col_species = col_species, genes = rownames(Y))
  class(res) <- "PhyloMash"
  return(res)
}

#' @title Call shifts at a Bayesian false discovery rate
#'
#' @description
#' Call the (gene, branch) pairs with the highest posterior probability of a
#' shift, as many as possible while keeping the expected proportion of false
#' calls (one minus the average posterior probability of the calls) below \code{level}.
#'
#' @param object a \code{PhyloMash} object, see \code{\link{phyloMash}}.
#' @param level the Bayesian false discovery rate. Default to 0.05.
#' @param by \code{"pair"} (the default) calls (gene, branch) pairs from \code{branch_prob};
#' \code{"gene"} calls genes with a shift on any branch from \code{shift_prob}.
#'
#' @return A logical matrix (genes by branches) or vector (genes).
#'
#' @export
#'
callShifts <- function(object, level = 0.05, by = c("pair", "gene")) {
  by <- match.arg(by)
  prob <- if (by == "pair") object$branch_prob else object$shift_prob
  p <- as.vector(prob)
  o <- order(-p)
  cum <- cumsum(1 - p[o]) / seq_along(o)
  k <- suppressWarnings(max(which(cum <= level)))
  sel <- rep(FALSE, length(p))
  if (is.finite(k)) sel[o[seq_len(k)]] <- TRUE
  if (by == "pair") {
    sel <- matrix(sel, nrow(prob), dimnames = dimnames(prob))
  } else names(sel) <- names(prob)
  sel
}
