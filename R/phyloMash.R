#' @title Phylogenetic shift detection on species means with mash
#'
#' @description
#' Detect branch-specific shifts in expression (or accessibility) across
#' species with replicates, by empirical Bayes on the species means. The
#' deviations of a gene's species means from their across-species mean are
#' modelled as sampling error (the standard errors of the species means, from
#' the replicates) plus phylogenetic drift plus, possibly, one shift shared by
#' all the species below one branch. The prior on drift and shifts is learned
#' from all genes. Two priors are available:
#' \itemize{
#' \item \code{prior = "factorized"} (the default): the prior weight of a
#' component is a product of independent factors, \eqn{u(d)\, q(s)\, [1 - \epsilon]}
#' without shift and \eqn{u(d)\, q(s)\, \epsilon_b\, h(v)} with a shift of
#' variance \eqn{v} on branch \eqn{b}, times \eqn{r(n)} for every component, where
#' \eqn{s} is the drift variance (per unit tree height), \eqn{d} the drift shape
#' (BM, or OU for the values of \code{rho_prime_grid}), \eqn{n} the variance of an
#' independent non-phylogenetic species-level noise (\code{noise_levels} times the
#' median squared standard error: outlier libraries, unmodelled batch), \eqn{h} a
#' distribution of shift sizes shared by all branches and \eqn{\epsilon_b} the
#' shift rate of branch \eqn{b} (one rate for all branches if
#' \code{shared_shift_rate}). The factors are estimated by EM. Because the
#' distribution of drift scales is the same for shifted and unshifted genes,
#' shifts cannot stand in for large drift (or the reverse), which, with few
#' species, the free weights of the mash prior allow at almost no cost in
#' likelihood; the Bayesian false discovery rate of \code{\link{callShifts}} is
#' then close to its nominal level, and as powerful as an oracle that knows the
#' population distributions of drift and shifts (package simulations with 4 and 7 species).
#' \item \code{prior = "mash"}: a \code{\link[mashr]{mash}} mixture of scaled
#' covariance patterns, each pattern and scale with its own weight: "drift"
#' patterns (BM, OU, drift plus noise, identity) and "shift" patterns (BM drift
#' plus a shift of variance \code{shift_ratios} times the drift variance, or a
#' pure shift), with the \code{drift_prior} penalty. A gene carries a shift on
#' branch \eqn{b} with the posterior weight of the shift patterns of \eqn{b} with
#' a ratio of at least \code{min_shift_ratio}.
#' }
#' See \code{\link{callShifts}} to call shifts at a given Bayesian false discovery rate.
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
#' @param noise_ratios mash prior: ratios of the species-specific non-phylogenetic variance
#' to the drift variance of the drift-plus-noise patterns.
#' @param identity mash prior: if \code{TRUE} (the default), add a pure species-specific variance pattern.
#' @param trend if \code{TRUE} (the default), the moderation of the within-species
#' variances allows for a mean-variance trend (\code{\link[limma]{squeezeVar}} with \code{covariate}).
#' @param prior \code{"factorized"} (the default) or \code{"mash"}, see Description.
#' @param shift_sd range of the shift standard deviations (log2 scale) of the factorized prior.
#' @param n_shift number of shift sizes of the factorized prior (log-spaced over \code{shift_sd}).
#' @param noise_levels factorized prior: variances of the non-phylogenetic species-level
#' noise, in multiples of the median squared standard error of the contrasts (the
#' level 0 is always included).
#' @param shared_shift_rate factorized prior: one shift rate for all branches (the
#' default; the rate is then estimated mostly from the branches where shifts and
#' drift are distinguishable, the short ones) or one rate per branch (more
#' flexible when shifts concentrate on some branches, slightly anti-conservative
#' in the package simulations).
#' @param max_iter maximum number of EM iterations of the factorized prior.
#' @param n_grid number of drift scales (factorized prior) or of scales of the mash grid.
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
#' @param fit_subset the number of (randomly chosen) genes used to estimate the
#' prior, which is then applied to all genes. \code{NULL}: all genes for the mash
#' prior, at most 5000 for the factorized prior.
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
#' \item{mash}{the \code{mash} fit (\code{prior = "mash"}).}
#' \item{prior_fit}{the estimated factors (\code{prior = "factorized"}): shift rate(s),
#' shift-size distribution, drift-scale distribution, drift-shape weights, noise distribution.}
#' \item{branch_prob}{genes by testable branches matrix of posterior probabilities of a shift.}
#' \item{shift_prob}{posterior probability of a shift on any branch.}
#' \item{top_branch}{the branch with the highest posterior probability.}
#' \item{pm_species}{genes by species matrix of posterior mean deviations from the across-species mean.}
#' \item{species_means, se_species}{the species means and their (moderated) standard errors.}
#' \item{pattern_weights}{the estimated prior weight of each pattern (summed over scales; mash prior),
#' or the shift rates (factorized prior).}
#' \item{branches}{the branch table (testable branches).}
#' }
#'
#' @seealso \code{\link{callShifts}}, \code{\link{branchIndicators}}, \code{\link[mashr]{mash}}
#'
#' @export
#'
phyloMash <- function(object, phy, col_species = NULL, model = c("BM", "OUfixedRoot"), prior = c("factorized", "mash"),
                      rho_prime_grid = c(0.1, 0.5), shift_sd = c(0.25, 5), n_shift = 8, noise_levels = c(0.1, 0.25, 0.5, 1, 2, 4),
                      shared_shift_rate = TRUE, max_iter = 1500,
                      shift_ratios = c(1, 4, 16, 64), min_shift_ratio = 4,
                      noise_ratios = c(0.25, 1, 4), identity = TRUE, trend = TRUE, n_grid = 15, drift_prior = 0.02,
                      fit_subset = NULL, covariates = NULL, patterns = c("tree", "canonical"), ...) {
  if (!requireNamespace("mashr", quietly = TRUE)) stop("Package 'mashr' is needed for function 'phyloMash'.")
  model <- match.arg(model)
  patterns <- match.arg(patterns)
  prior <- if (patterns == "canonical") "mash" else match.arg(prior)
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
  drift_vcv <- list(BM = ape::vcv(tr_u))
  if (model == "OUfixedRoot") {
    ab <- getBoundsSelectionStrength(tr_u)
    for (rp in rho_prime_grid) drift_vcv[[paste0("OU", rp)]] <- ape::vcv(rescale_tree(phylolm::transf.branch.lengths(tr_u, "OUfixedRoot", parameters = list(alpha = rho_prime_inv(rp, 1, ab)))$tree))
  }
  if (prior == "factorized") {
    ## drift shapes (BM, OU), an independent non-phylogenetic noise factor (iid species-level variance, e.g. from
    ## outlier libraries or unmodelled batch, on its own grid), and the shift factor
    shapes <- lapply(drift_vcv, function(Cd) { x <- L %*% Cd %*% t(L); (x + t(x)) / 2 })
    noise_shape <- L %*% t(L)
    dirs <- lapply(seq_along(labels), function(b) as.vector(L %*% Zs[, b]))
    ## drift variance per unit tree height: from a hundredth of the smallest standard error to four times the largest squared deviation
    s_hi <- 4 * stats::quantile(rowMeans(Bhat^2), 0.999) / mean(diag(shapes$BM)); s_lo <- (min(Shat) / 10)^2 / mean(diag(shapes$BM))
    s_grid <- c(0, exp(seq(log(s_lo), log(s_hi), length.out = n_grid)))
    v_grid <- exp(seq(log(shift_sd[1]^2), log(shift_sd[2]^2), length.out = n_shift))
    n_grid_noise <- c(0, noise_levels * stats::median(data$Shat^2))
    data_fit <- data
    n_fit <- if (is.null(fit_subset)) min(G, 5000) else min(G, fit_subset)
    if (n_fit < G) {
      sub <- sort(sample.int(G, n_fit))
      data_fit <- mashr::mash_set_data(data$Bhat[sub, , drop = FALSE], data$Shat[sub, , drop = FALSE], V = data$V)
    }
    if (n_fit < 1000) warning("The factorized prior is fitted on ", n_fit, " genes only: its shift rate may be poorly estimated (it has more components than genes).")
    fit <- fit_factorized_prior(data_fit, shapes, noise_shape, dirs, s_grid, n_grid_noise, v_grid, shared = shared_shift_rate, max_iter = max_iter)
    post <- factorized_posterior(data, fit)
    branch_prob <- post$branch_prob; dimnames(branch_prob) <- list(stats$genes, labels)
    pm_species <- cbind(post$pm, -rowSums(post$pm)); colnames(pm_species) <- species
    e <- fit$e; names(e) <- c("no_shift", labels)
    res <- list(branch_prob = branch_prob, shift_prob = rowSums(branch_prob),
                top_branch = labels[max.col(branch_prob, ties.method = "first")],
                pm_species = pm_species, species_means = M, se_species = SE, dropout = stats$dropout[, ord, drop = FALSE],
                pattern_weights = e,
                prior_fit = list(shift_rate = e, shift_size = data.frame(sd = sqrt(v_grid), weight = fit$h),
                                 drift_scale = data.frame(variance = s_grid, weight = fit$q), drift_shape = setNames(fit$u, names(shapes)),
                                 noise = data.frame(variance = n_grid_noise, weight = fit$r),
                                 loglik = fit$loglik, iterations = fit$iter),
                branches = branches, L = L, df = df, data = data, prior = prior,
                model = model, patterns = patterns, phy = phy, stats_method = stats$method, genes = stats$genes)
    class(res) <- "PhyloMash"
    return(res)
  }
  toc <- function(U) { Uc <- L %*% U %*% t(L); Uc <- (Uc + t(Uc)) / 2; Uc / max(diag(Uc)) }
  Ulist <- list(); tag <- character()
  if (patterns == "tree") {
    drifts <- drift_vcv
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
              data = data, tag = tag, prior = prior,
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

## Factorized prior: component = (drift shape d, drift scale s) x {no shift, shift on branch b with variance v};
## pi = u[d] q[s] (e0 | e[b] h[v]). Maximum likelihood by EM (product-form mixture: each factor is the marginal of the
## responsibilities), accelerated by SQUAREM (Varadhan & Roland 2008), with a fallback to plain EM steps.
fit_factorized_prior <- function(data, shapes, noise_shape, dirs, s_grid, n_grid, v_grid, shared = TRUE, max_iter = 1500, tol = 1e-10) {
  nd <- length(shapes); ns <- length(s_grid); nn <- length(n_grid); nb <- length(dirs); nv <- length(v_grid)
  comp <- expand.grid(d = seq_len(nd), s = seq_len(ns), n = seq_len(nn), j = 0:(nb * nv))
  comp$b <- ifelse(comp$j == 0, 0L, (comp$j - 1L) %/% nv + 1L); comp$v <- ifelse(comp$j == 0, 0L, (comp$j - 1L) %% nv + 1L)
  U <- lapply(seq_len(nrow(comp)), function(i) { x <- s_grid[comp$s[i]] * shapes[[comp$d[i]]] + n_grid[comp$n[i]] * noise_shape
    if (comp$b[i] > 0) x <- x + v_grid[comp$v[i]] * tcrossprod(dirs[[comp$b[i]]]); (x + t(x)) / 2 })
  lik <- exp(mashr:::calc_relative_lik_matrix(data, U)$loglik_matrix); G <- nrow(lik)
  blocks <- c(rep("u", nd), rep("q", ns), rep("r", nn), rep("e", nb + 1), rep("h", nv))
  unpack <- function(th) list(u = th[blocks == "u"], q = th[blocks == "q"], r = th[blocks == "r"], e = th[blocks == "e"], h = th[blocks == "h"])
  pi_of <- function(p) p$u[comp$d] * p$q[comp$s] * p$r[comp$n] * p$e[comp$b + 1L] * ifelse(comp$b == 0, 1, p$h[pmax(comp$v, 1L)])
  loglik <- function(th) sum(log(as.vector(lik %*% pi_of(unpack(th)))))
  em <- function(th) {
    p <- unpack(th); pic <- pi_of(p); f <- as.vector(lik %*% pic)
    cs <- pic * as.vector(crossprod(lik, 1 / f)) / G
    u <- as.vector(tapply(cs, comp$d, sum)); q <- as.vector(tapply(cs, comp$s, sum)); r <- as.vector(tapply(cs, comp$n, sum))
    e <- as.vector(tapply(cs, comp$b, sum)); if (shared) e <- c(e[1], rep(sum(e[-1]) / nb, nb))
    h <- as.vector(tapply(cs[comp$b > 0], comp$v[comp$b > 0], sum)); h <- h / sum(h)
    c(u, q, r, e, h) }
  th <- c(rep(1 / nd, nd), rep(1 / ns, ns), rep(1 / nn, nn), c(0.9, rep(0.1 / nb, nb)), rep(1 / nv, nv))
  ## SQUAREM on the log scale (weights stay positive), each block renormalized; weights floored to keep the logs finite
  lg <- function(x) log(pmax(x, 1e-300))
  renorm <- function(x) { for (bk in unique(blocks)) x[blocks == bk] <- x[blocks == bk] / sum(x[blocks == bk]); x }
  ll <- loglik(th); n_em <- 0
  while (n_em < max_iter) {
    th1 <- em(th); th2 <- em(th1); n_em <- n_em + 2
    r <- lg(th1) - lg(th); v <- lg(th2) - lg(th1) - r
    ll2 <- loglik(th2)
    if (sum(v^2) < 1e-30) { th_new <- th2; ll_new <- ll2 } else {
      alpha <- min(-1, -sqrt(sum(r^2) / sum(v^2)))
      thp <- renorm(exp(pmin(lg(th) - 2 * alpha * r + alpha^2 * v, 700)))
      thp <- em(thp); n_em <- n_em + 1; ll_p <- loglik(thp)
      if (is.finite(ll_p) && ll_p >= ll2) { th_new <- thp; ll_new <- ll_p } else { th_new <- th2; ll_new <- ll2 }
    }
    converged <- abs(ll_new - ll) < tol * abs(ll) && max(abs(th_new - th)) < 1e-6
    th <- th_new; ll <- ll_new
    if (converged) break
  }
  p <- unpack(th)
  list(comp = comp, U = U, pi = pi_of(p), u = p$u, q = p$q, r = p$r, e = p$e, h = p$h, loglik = ll, iter = n_em)
}

## Posterior shift probabilities and posterior mean contrasts for all genes, in chunks
factorized_posterior <- function(data, fit, chunk = NULL) {
  G <- nrow(data$Bhat); nb <- max(fit$comp$b)
  use <- which(fit$pi > 1e-12); U <- fit$U[use]; pi <- fit$pi[use]; bcomp <- fit$comp$b[use]
  if (is.null(chunk)) chunk <- max(200, floor(5e6 / length(use)))
  bp <- matrix(0, G, nb); pm <- matrix(0, G, ncol(data$Bhat))
  for (st in seq(1, G, by = chunk)) {
    idx <- st:min(G, st + chunk - 1)
    dc <- mashr::mash_set_data(data$Bhat[idx, , drop = FALSE], data$Shat[idx, , drop = FALSE], V = data$V)
    lik <- exp(mashr:::calc_relative_lik_matrix(dc, U)$loglik_matrix)
    w <- sweep(lik, 2, pi, "*"); w <- w / rowSums(w)
    bp[idx, ] <- sapply(seq_len(nb), function(b) rowSums(w[, bcomp == b, drop = FALSE]))
    pm[idx, ] <- mashr:::compute_posterior_matrices(dc, U, w)$PosteriorMean
  }
  list(branch_prob = bp, pm = pm)
}
