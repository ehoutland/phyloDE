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
#' @param object a matrix of normalized log2 expression (genes by samples), or a
#' \code{SpeciesStats} object from \code{\link{speciesStats}} (e.g. with
#' \code{method = "edgeR"}, for read-depth-aware, species-specific standard errors).
#' A matrix is passed to \code{\link{speciesStats}} with \code{method = "limma"}.
#' @param phy the species tree (one tip per species), dated and ultrametric.
#' @param model the phylogenetic model of the drift: \code{"BM"} (the default) or
#' \code{"OUfixedRoot"}.
#' @param covariates optional sample-level covariates, see \code{\link{speciesStats}}
#' (ignored if \code{object} is a \code{SpeciesStats} object).
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
#' @param drift_prior pseudo-observations explained by the no-shift patterns only
#' (null, drift, drift plus noise, identity, and shift patterns below
#' \code{min_shift_ratio}), as a fraction of the number of genes used to fit the
#' mixture weights; see \code{\link{fit_mash_weights}}. With few species the
#' drift and drift-plus-shift patterns explain most genes equally well and the
#' fitted weight of the shift patterns, hence the posterior shift probabilities
#' of unshifted genes, is inflated; the pseudo-observations resolve this tie in
#' favour of no shift. The default \code{0.02} keeps the realized false discovery
#' rate of \code{\link{callShifts}} at or below its nominal level up to 0.05 in
#' simulations with 4 to 7 species, for a loss of about 0.05 in power; \code{0}
#' is the plain \code{mash} fit, which is anti-conservative with few species
#' (realized 0.07-0.10 at nominal 0.05 on branches, 0.13-0.14 on genes with 4 species).
#' @param fit_subset if not \code{NULL}, the number of (randomly chosen) genes used to
#' estimate the mixture weights, which are then applied to all genes
#' (faster for large data sets).
#' @param patterns \code{"tree"} (the default): the drift and shift patterns described
#' above; \code{"canonical"}: the data-agnostic patterns of \code{\link[mashr]{cov_canonical}}
#' only (null, singletons, equal effects, simple heterogeneity, identity), as in
#' \code{dreamlet::run_mash}, with the null correlation estimated by
#' \code{\link[mashr]{mash_estimate_corr_em}}; a singleton pattern on species s counts
#' as a shift on the terminal branch of s, and internal branches get no shift
#' probability. Useful as a reference for what the tree-derived patterns add.
#' @param ... further arguments passed to \code{\link{fit_mash_weights}} (e.g. \code{control} for \code{mixsqp}).
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
                      noise_ratios = c(0.25, 1, 4), identity = TRUE, trend = TRUE, n_grid = 15, drift_prior = 0.02,
                      fit_subset = NULL, covariates = NULL, patterns = c("tree", "canonical"), ...) {
  if (!requireNamespace("mashr", quietly = TRUE)) stop("Package 'mashr' is needed for function 'phyloMash'.")
  model <- match.arg(model)
  patterns <- match.arg(patterns)
  if (!ape::is.ultrametric(phy)) stop("The tree must be ultrametric.")
  if (inherits(object, "SpeciesStats")) {
    stats <- object
  } else {
    y <- check_expression_matrix(object)
    if (is.null(col_species)) col_species <- parse_species(phy, colnames(y$exprs))
    stats <- speciesStats(y$exprs, col_species = col_species, method = "limma", covariates = covariates, trend = trend)
  }
  if (!all(stats$species %in% phy$tip.label)) stop("Some species of the data are not in the tree.")
  if (!all(phy$tip.label %in% stats$species)) stop("All the species of the tree must have samples.")
  ord <- match(phy$tip.label, stats$species)
  M <- stats$mean[, ord, drop = FALSE]; SE <- stats$se[, ord, drop = FALSE]
  species <- phy$tip.label
  S <- length(species); G <- nrow(M)
  df <- stats$df
  if (length(df) > 1) df <- matrix(df, G, S - 1)

  ## contrasts: deviations from the unweighted across-species mean, with gene- and species-specific errors
  L <- mashr::contrast_matrix(S, "mean", species)
  Bhat <- M %*% t(L)
  Shat <- sqrt(SE^2 %*% t(L^2))
  ## null correlation of the contrasts: depends on the gene through SE; use the average measurement variances
  CL <- L %*% diag(colMeans(SE^2) / mean(SE^2)) %*% t(L)
  V <- stats::cov2cor(CL)
  data <- mashr::mash_set_data(Bhat, Shat, df = df, V = V)

  ## covariance patterns in species space, mapped to the contrast space
  bi <- branchIndicators(phy, species)
  keep <- bi$branches$testable
  Zs <- bi$Z[, keep, drop = FALSE]
  branches <- bi$branches[keep, ]
  labels <- branches$label
  tr_u <- rescale_tree(phy)
  toc <- function(U) { Uc <- L %*% U %*% t(L); Uc <- (Uc + t(Uc)) / 2; Uc / max(diag(Uc)) }
  Ulist <- list(); tag <- character()
  if (patterns == "tree") {
    drifts <- list(BM = ape::vcv(tr_u))
    if (model == "OUfixedRoot") {
      ab <- getBoundsSelectionStrength(tr_u)
      for (rp in rho_prime_grid) {
        al <- rho_prime_inv(rp, 1, ab)
        drifts[[paste0("OU", rp)]] <- ape::vcv(rescale_tree(phylolm::transf.branch.lengths(tr_u, "OUfixedRoot", parameters = list(alpha = al))$tree))
      }
    }
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
  } else {
    ## dreamlet-style: canonical patterns in the contrast space; singleton on contrast "s - mean" ~ shift on species s
    Ulist <- mashr::cov_canonical(data)
    tag <- setNames(rep("none", length(Ulist)), names(Ulist))
    for (s_i in seq_len(S - 1)) {
      nm <- paste0("singletons_", s_i)
      if (nm %in% names(Ulist)) tag[nm] <- species[s_i]
    }
    ## the S-th species is -sum of the contrasts: no singleton; use the equal-effects-with-opposite-sign pattern as its shift
    Ulist[[paste0("singleton_", species[S])]] <- toc(tcrossprod(as.numeric(seq_len(S) == S)))
    tag[paste0("singleton_", species[S])] <- species[S]
  }

  ## mash fit: mixture weights on a (sub)set of genes, with the no-shift pseudo-observations, then posteriors for all genes
  g_auto <- mashr:::autoselect_grid(data, sqrt(2))
  grid <- exp(seq(log(min(g_auto)), log(max(g_auto)), length.out = n_grid))
  if (patterns == "canonical") {
    ## estimate the null correlation from the data (dreamlet::run_mash)
    Vem <- mashr::mash_estimate_corr_em(data, Ulist, details = TRUE)
    data <- mashr::mash_update_data(data, V = Vem$V)
  }
  data_fit <- data
  if (!is.null(fit_subset) && fit_subset < G) {
    sub <- sort(sample.int(G, fit_subset))
    df_sub <- if (is.matrix(df)) df[sub, , drop = FALSE] else df
    data_fit <- mashr::mash_set_data(Bhat[sub, , drop = FALSE], Shat[sub, , drop = FALSE], df = df_sub, V = data$V)
  }
  g <- fit_mash_weights(data_fit, Ulist, grid, no_shift = c("null", names(tag)[tag == "none"]), pseudo = drift_prior, ...)
  m <- mashr::mash(data, g = g, fixg = TRUE, outputlevel = 2, verbose = FALSE)

  ## posterior probabilities of the shift patterns
  branch_prob <- branch_prob_from_weights(m$posterior_weights, tag, labels)
  dimnames(branch_prob) <- list(stats$genes, labels)
  pi <- m$fitted_g$pi
  pattern_weights <- tapply(pi, sub("\\.[0-9]+$", "", names(pi)), sum)
  pmc <- ashr::get_pm(m)
  pm_species <- cbind(pmc, -rowSums(pmc)); colnames(pm_species) <- species

  res <- list(mash = m, branch_prob = branch_prob, shift_prob = rowSums(branch_prob),
              top_branch = labels[max.col(branch_prob, ties.method = "first")],
              pm_species = pm_species, species_means = M, se_species = SE, dropout = stats$dropout[, ord, drop = FALSE],
              pattern_weights = pattern_weights, branches = branches, L = L, df = df, drift_prior = drift_prior,
              data = data, tag = tag,
              model = model, patterns = patterns, phy = phy, stats_method = stats$method, genes = stats$genes)
  class(res) <- "PhyloMash"
  return(res)
}

#' @title Fit the mash mixture weights with pseudo-observations on a group of components
#'
#' @description
#' The weights of a \code{\link[mashr]{mash}} mixture are the maximizer of the
#' marginal likelihood of the genes, with a Dirichlet prior (\code{mash}'s
#' default: \code{nullweight} pseudo-observations on the null point mass). Here
#' \code{pseudo} times the number of genes extra pseudo-observations are added
#' that can be explained by any of the \code{no_shift} components (null, drift,
#' ...) but by none of the others: a Dirichlet prior on the split between the
#' "no shift" group and the shift components. When the likelihood is flat
#' between the two (few species: drift and drift-plus-shift patterns explain
#' most genes equally well), the penalty resolves the tie in favour of no
#' shift; genes carrying a real shift still pull their weight. The weights are
#' then fixed for the posterior computation (\code{mash(..., fixg = TRUE)}).
#'
#' @param data a mash data object (the genes used to fit the weights).
#' @param Ulist,grid the covariance patterns and scales (as in \code{mash}).
#' @param no_shift names of the components counting as "no shift" (pattern
#' names of \code{Ulist}, and \code{"null"} for the point mass).
#' @param pseudo pseudo-observations on the no-shift group, as a fraction of the
#' number of genes in \code{data}.
#' @param nullweight prior weight of the null point mass, as in \code{mash} (default 10).
#' @param control passed to \code{\link[mixsqp]{mixsqp}}.
#' @param ... ignored (so that \code{mash} arguments can be passed through).
#'
#' @return A fitted \code{g} (list \code{pi}, \code{Ulist}, \code{grid}, \code{usepointmass}).
#'
#' @keywords internal
#'
fit_mash_weights <- function(data, Ulist, grid, no_shift = "null", pseudo = 0, nullweight = 10, control = list(), ...) {
  xU <- mashr:::expand_cov(Ulist, grid, usepointmass = TRUE)
  L <- exp(mashr:::calc_relative_lik_matrix(data, xU)$loglik_matrix)
  K <- ncol(L); J <- nrow(L)
  base <- c("null", sub("\\.[0-9]+$", "", names(xU)[-1]))
  A <- rbind(as.numeric(seq_len(K) == 1), L); w <- c(nullweight - 1, rep(1, J))
  if (pseudo > 0) { A <- rbind(A, as.numeric(base %in% no_shift)); w <- c(w, pseudo * J) }
  keep <- w > 0
  ctrl <- utils::modifyList(list(verbose = FALSE, eps = 1e-6, numiter.em = 20), control)
  out <- mixsqp::mixsqp(A[keep, , drop = FALSE], w[keep], control = ctrl)
  pi <- out$x / sum(out$x); names(pi) <- names(xU)
  list(pi = pi, Ulist = Ulist, grid = grid, usepointmass = TRUE)
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
  prob <- if (by == "pair") object$branch_prob else pmin(object$shift_prob, 1)
  sel <- bfdr_select(as.vector(prob), level)
  if (by == "pair") sel <- matrix(sel, nrow(prob), dimnames = dimnames(prob)) else names(sel) <- names(prob)
  sel
}

## Bayesian FDR selection: largest set of top scores whose average (1 - p) is <= level
bfdr_select <- function(p, level) {
  o <- order(-p); cum <- cumsum(1 - p[o]) / seq_along(o)
  k <- suppressWarnings(max(which(cum <= level)))
  sel <- rep(FALSE, length(p)); if (is.finite(k)) sel[o[seq_len(k)]] <- TRUE
  sel
}

## Shift probabilities per branch from mash posterior weights (columns named by component)
branch_prob_from_weights <- function(w, tag, labels) {
  base <- sub("\\.[0-9]+$", "", colnames(w))
  tag_w <- tag[base]; tag_w[base == "null"] <- "none"
  bp <- sapply(labels, function(b) rowSums(w[, which(tag_w == b), drop = FALSE]))
  if (is.null(dim(bp))) bp <- matrix(bp, nrow = nrow(w))
  bp
}
