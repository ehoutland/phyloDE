#' @title Precision weights from counts (voom-style)
#'
#' @description
#' Observation-level precision weights for count data, following
#' \code{\link[limma]{voom}} and \code{voomWithDreamWeights} of the
#' \code{variancePartition} package: a lowess trend of the square-root
#' residual standard deviation on the mean log-count is fitted across genes,
#' and the weight of each observation is the inverse of the squared trend
#' evaluated at its own fitted log-count. A sample with few reads for a gene
#' (small library, or low expression in its species) therefore gets a small
#' weight, a sample with many reads a large one.
#'
#' @param counts genes by samples matrix of counts.
#' @param col_species the species of each sample (the trend uses species means as fitted values).
#' @param lib.size library sizes (default: column sums times \code{norm.factors}).
#' @param norm.factors normalization factors (e.g. TMM); default 1.
#' @param span lowess span. Default 0.5.
#' @param prior.count added to the counts before the log2 CPM. Default 0.5.
#'
#' @return A list with \code{E} (log2 CPM), \code{weights} (genes by samples) and
#' the trend (\code{trend}, a function).
#'
#' @export
#'
voomPhyloWeights <- function(counts, col_species, lib.size = NULL, norm.factors = NULL, span = 0.5, prior.count = 0.5) {
  counts <- as.matrix(counts)
  n <- ncol(counts); G <- nrow(counts)
  if (is.null(norm.factors)) norm.factors <- rep(1, n)
  if (is.null(lib.size)) lib.size <- colSums(counts) * norm.factors
  E <- t(log2(t(counts + prior.count) / (lib.size + 1) * 1e6))
  sp <- factor(col_species)
  X <- stats::model.matrix(~ 0 + sp)
  ## species-mean fit (fixed effects only, as voom uses the design fit)
  fit <- limma::lmFit(E, X)
  fitted_logcpm <- fit$coefficients %*% t(X)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e6)
  sy <- sqrt(fit$sigma)
  ok <- is.finite(sx) & is.finite(sy)
  l <- stats::lowess(sx[ok], sy[ok], f = span)
  f <- stats::approxfun(l, rule = 2)
  fitted_count <- t(t(2^fitted_logcpm) * (lib.size + 1)) * 1e-6
  fitted_logcount <- log2(fitted_count)
  w <- 1 / f(fitted_logcount)^4
  dim(w) <- dim(counts); dimnames(w) <- dimnames(counts)
  dimnames(E) <- dimnames(counts)
  list(E = E, weights = w, trend = f, lib.size = lib.size)
}

#' @title Phylogenetic linear model with precision weights (weighted phylogenetic GLS)
#'
#' @description
#' For every gene, fit \eqn{y = X\beta + u + e}, with \eqn{u \sim N(0, \sigma^2_u C(\theta))}
#' the phylogenetic (species-level, tree-correlated) random effect and
#' \eqn{e \sim N(0, \sigma^2_e \, \mathrm{diag}(1 / w))} heteroscedastic
#' sample-level noise with known precision weights \eqn{w} (from
#' \code{\link{voomPhyloWeights}}). This is the model of \code{dream}
#' (\code{variancePartition}) with the \code{(1 | species)} random effect
#' replaced by the phylogenetic covariance, and the gene-wise variance
#' components handled on the empirical Bayes grid of \code{\link{phyloEB}}
#' instead of a per-gene REML optimization: \eqn{\lambda = \sigma^2_u / (\sigma^2_u + \sigma^2_e / \bar{w})}
#' and \eqn{\rho'} (OU) are gridded, the scale is profiled, the prior over the
#' grid is estimated across genes and the estimates are averaged over each
#' gene's posterior (law of total variance for the standard errors).
#'
#' Tests: the last coefficient of \code{design} (if \code{design} has more than
#' an intercept), and/or a shift on every testable branch of the species tree
#' (indicator of the descendant samples added to \code{design}), as in
#' \code{\link{phyloShiftFit}}.
#'
#' @param counts genes by samples matrix of counts (or a list with \code{E} and
#' \code{weights} as returned by \code{\link{voomPhyloWeights}}).
#' @param phy the species tree (one tip per species), ultrametric.
#' @param col_species the species of each sample.
#' @param design samples by coefficients design matrix with an intercept; the
#' last coefficient is tested. Default: intercept only.
#' @param shifts if \code{TRUE} (the default), test a shift on every testable branch.
#' @param model \code{"BM"} (the default) or \code{"OUfixedRoot"}.
#' @param n_lambda,n_rho,n_eta grid sizes for \eqn{\lambda}, \eqn{\rho'} and \eqn{\eta}.
#' @param species_variance if \code{TRUE} (the default), the model includes a third
#' variance component, a species-level non-phylogenetic random effect
#' \eqn{v \sim N(0, \sigma^2_v I_{species})}, with \eqn{\eta = \sigma^2_v / (\sigma^2_u + \sigma^2_v)}
#' the non-phylogenetic share of the species-level variance on a grid
#' (\code{n_eta} values). Depth- or batch-related species effects that do not
#' follow the tree (e.g. systematic sequencing-depth differences between species
#' at low counts) are absorbed by this term instead of being read as shifts.
#' @param eta_null_weight the prior over the grid is estimated with a null-biased
#' penalty on \eqn{\eta} (pseudo-observations explained only by \eqn{\eta = 0}
#' components, as the \code{nullbiased} prior of \code{ashr}), with this weight
#' relative to the number of genes. Default 0 (unpenalised): in simulations, the
#' penalty only reduced the protection against species-level nuisance variance.
#' @param weights optional genes by samples matrix of precision weights, used
#' instead of \code{\link{voomPhyloWeights}}. \code{NULL} with a count matrix
#' computes voom-style weights; a matrix of ones gives the unweighted fit.
#' @param norm.factors normalization factors for \code{\link{voomPhyloWeights}}.
#' @param squeeze moderate the residual variances across genes with
#' \code{\link[limma]{squeezeVar}} (default \code{TRUE}).
#' @param n_profiles genes are grouped by k-means on their (log) normalized
#' inverse-weight profiles into at most this many groups, and the genes of a
#' group share the Cholesky factors of their group centroid (default 200;
#' \code{Inf} for exact gene-wise factors).
#' @param prior_smoothing,prior_pseudocount see \code{\link{phyloEB}}.
#'
#' @return An object of class \code{PhyloShiftFit} (see \code{\link{phyloShiftFit}}),
#' with the design coefficient as the first column when tested. Also
#' \code{weights} and \code{E}.
#'
#' @seealso \code{\link{voomPhyloWeights}}, \code{\link{phyloShiftFit}}, \code{\link{phyloShiftGLS}}
#'
#' @export
#'
phyloWeightedFit <- function(counts, phy, col_species, design = NULL, shifts = TRUE,
                             model = c("BM", "OUfixedRoot"), n_lambda = 15, n_rho = 7, n_eta = 5, species_variance = TRUE, eta_null_weight = 0,
                             weights = NULL, norm.factors = NULL, squeeze = TRUE, n_profiles = 200,
                             prior_smoothing = c(0.5, 1, 2), prior_pseudocount = 1e-8) {
  model <- match.arg(model)
  if (is.list(counts) && !is.null(counts$E)) { E <- counts$E; W <- counts$weights }
  else {
    v <- voomPhyloWeights(counts, col_species, norm.factors = norm.factors)
    E <- v$E; W <- if (is.null(weights)) v$weights else weights
  }
  if (!is.null(weights)) W <- weights
  n <- ncol(E); G <- nrow(E)
  y <- check_expression_matrix(E)
  phy_rep <- check_tree(phy, y, col_species)
  Y <- checkParamMatrix(y$exprs, "expression matrix", phy_rep)
  W <- W[, colnames(Y), drop = FALSE]
  if (is.null(design)) design <- matrix(1, n, 1, dimnames = list(colnames(E), "(Intercept)"))
  design <- check_design_matrix(design, y, phy_rep)
  sample_species <- col_species[match(colnames(Y), colnames(E))]

  tests <- list()
  if (ncol(design) > 1) tests[[colnames(design)[ncol(design)]]] <- list(X0 = design[, -ncol(design), drop = FALSE], z = design[, ncol(design)])
  branches <- NULL; Z <- NULL
  if (shifts) {
    bi <- branchIndicators(phy, sample_species)
    keep <- bi$branches$testable
    estimable <- sapply(seq_len(nrow(bi$branches)), function(b) qr(cbind(design, bi$Z[, b]))$rank > ncol(design))
    keep <- keep & estimable
    for (b in which(keep)) tests[[bi$branches$label[b]]] <- list(X0 = design, z = bi$Z[, b])
    branches <- bi$branches[keep, ]; Z <- bi$Z[, keep, drop = FALSE]; rownames(Z) <- colnames(Y)
  }

  ## grid
  tr_u <- rescale_tree(phy_rep)
  lambda_grid <- c(0.001, stats::plogis(seq(-5, 5, length.out = n_lambda - 2)), 0.999)
  eta_grid <- if (species_variance) c(0, stats::plogis(seq(-3, 3, length.out = n_eta - 1))) else 0
  if (model == "OUfixedRoot") {
    ab <- getBoundsSelectionStrength(tr_u)
    rb <- rho_prime(rev(ab), 1); rb <- c(max(rb[1], 1e-3), min(rb[2], 1 - 1e-3))
    rho_grid <- stats::plogis(seq(stats::qlogis(rb[1]) + 0.05, stats::qlogis(rb[2]) - 0.05, length.out = n_rho))
    grid <- expand.grid(lambda = lambda_grid, rho_prime = rho_grid, eta = eta_grid)
    grid$alpha <- sapply(grid$rho_prime, function(r) rho_prime_inv(r, 1, ab))
  } else grid <- expand.grid(lambda = lambda_grid, rho_prime = 1, eta = eta_grid, alpha = 0)
  grid$rho <- 1 - grid$rho_prime
  K <- nrow(grid)
  ## species-block identity on the replicate tree (samples of the same species share a block)
  I_sp <- outer(sample_species, sample_species, "==") * 1
  Cs <- lapply(seq_len(K), function(k) (1 - grid$eta[k]) * ape::vcv(eb_transform_tree(tr_u, model, 1, grid$alpha[k])) + grid$eta[k] * I_sp)
  ## normalized inverse weights (mean 1 per gene)
  Dn <- (1 / W) / rowMeans(1 / W)
  cst <- log(2 * pi) + 1

  ## genes with similar weight profiles share a Cholesky: cluster the normalized inverse-weight
  ## profiles (log scale) into at most `n_profiles` groups and use the group centroid
  n_profiles <- min(n_profiles, G, nrow(unique(round(log(Dn), 6))))
  if (is.finite(n_profiles) && G > n_profiles) {
    km <- stats::kmeans(log(Dn), centers = n_profiles, iter.max = 30, nstart = 1)
    prof <- exp(km$centers); grp <- km$cluster
  } else { prof <- Dn; grp <- seq_len(G) }
  ngrp <- nrow(prof)

  run_test <- function(X0, z) {
    p <- ncol(X0) + 1
    LL <- est <- s2 <- nz <- matrix(NA_real_, G, K)
    for (k in seq_len(K)) {
      Ck <- Cs[[k]]; lam <- grid$lambda[k]
      for (j in seq_len(ngrp)) {
        genes <- which(grp == j)
        Sig <- lam * Ck + (1 - lam) * diag(prof[j, ], n)
        R <- chol(Sig)
        X0t <- backsolve(R, X0, transpose = TRUE); zt <- backsolve(R, z, transpose = TRUE)
        q0 <- qr(X0t); zp <- qr.resid(q0, zt); nzk <- sum(zp^2)
        Yt <- backsolve(R, t(Y[genes, , drop = FALSE]), transpose = TRUE)
        Yp <- qr.resid(q0, Yt)
        b <- as.vector(crossprod(Yp, zp)) / nzk
        RSS <- colSums(Yp^2) - b^2 * nzk
        est[genes, k] <- b; nz[genes, k] <- nzk; s2[genes, k] <- RSS / (n - p)
        LL[genes, k] <- -0.5 * ((n - p) * (cst + log(RSS / (n - p))) + 2 * sum(log(diag(R))) + 2 * sum(log(abs(diag(qr.R(qr(cbind(X0t, zt))))))))
      }
    }
    pp <- eb_fit_prior(LL, grid, prior_smoothing, prior_pseudocount, null_weight = if (species_variance) eta_null_weight else 0)
    Wp <- pp$posterior
    if (squeeze) {
      sq <- limma::squeezeVar(rowSums(Wp * s2), df = n - p)
      df.prior <- if (is.finite(sq$df.prior)) sq$df.prior else 0
      s2 <- (df.prior * sq$var.prior + (n - p) * s2) / (df.prior + n - p)
    } else df.prior <- 0
    df <- n - p + df.prior
    se2 <- s2 / nz
    Bhat <- rowSums(Wp * est)
    Shat <- sqrt(rowSums(Wp * (se2 + (est - Bhat)^2)))
    pv <- rowSums(Wp * 2 * stats::pt(-abs(est) / sqrt(se2), df = df))
    list(Bhat = Bhat, Shat = Shat, p = pv, df = df, lambda = as.vector(Wp %*% grid$lambda), rho = as.vector(Wp %*% grid$rho), eta = as.vector(Wp %*% grid$eta), prior = pp$prior)
  }
  out <- lapply(tests, function(tt) run_test(tt$X0, tt$z))
  nm <- names(tests)
  mk <- function(f) { m <- sapply(out, `[[`, f); if (is.null(dim(m))) m <- matrix(m, G); dimnames(m) <- list(rownames(Y), nm); m }
  Bhat <- mk("Bhat"); Shat <- mk("Shat")
  res <- list(Bhat = Bhat, Shat = Shat, t = Bhat / Shat, p.value = mk("p"), df = out[[1]]$df,
              lambda = mk("lambda"), rho = mk("rho"), eta = mk("eta"), nuisance = "eb_weighted", branches = branches, Z = Z,
              priors = sapply(out, `[[`, "prior"), grid = grid, model = model, phy = phy, phy_rep = phy_rep,
              weights = W, E = Y, design = design, col_species = sample_species, genes = rownames(Y))
  class(res) <- "PhyloShiftFit"
  res
}

#' @title Weighted species-level statistics
#'
#' @description
#' Species means of the log2 CPM weighted by voom-style precision weights
#' (\code{\link{voomPhyloWeights}}), with standard errors from the weights and a
#' moderated residual scale: the heteroscedastic counterpart of
#' \code{speciesStats(method = "limma")}, for \code{\link{phyloMash}}.
#'
#' @inheritParams voomPhyloWeights
#' @param trend moderate the residual scale with a mean trend (default \code{TRUE}).
#'
#' @return A \code{SpeciesStats} object (see \code{\link{speciesStats}}), \code{method = "weighted"}.
#'
#' @export
#'
speciesStatsWeighted <- function(counts, col_species, norm.factors = NULL, trend = TRUE, span = 0.5, prior.count = 0.5) {
  v <- voomPhyloWeights(counts, col_species, norm.factors = norm.factors, span = span, prior.count = prior.count)
  E <- v$E; W <- v$weights
  sp <- factor(col_species); S <- nlevels(sp); n <- ncol(E); G <- nrow(E)
  Xs <- stats::model.matrix(~ 0 + sp)
  M <- (E * W) %*% Xs / (W %*% Xs)                      # weighted species means
  colnames(M) <- levels(sp)
  resid <- E - M[, as.integer(sp), drop = FALSE]
  s2w <- rowSums(W * resid^2) / (n - S)                   # weighted residual scale (phi in e ~ N(0, phi / w))
  sq <- if (trend) limma::squeezeVar(s2w, df = n - S, covariate = rowMeans(E)) else limma::squeezeVar(s2w, df = n - S)
  s2post <- sq$var.post
  df <- n - S + if (is.finite(sq$df.prior)) sq$df.prior else 0
  se <- sqrt(s2post / (W %*% Xs))                          # var(weighted mean) = phi / sum(w)
  colnames(se) <- levels(sp)
  dropout <- (counts %*% Xs) == 0; dimnames(dropout) <- dimnames(M)
  res <- list(mean = M, se = se, df = df, dropout = dropout, n_rep = setNames(as.vector(table(sp)), levels(sp)),
              method = "weighted", species = levels(sp), genes = rownames(E), weights = W, E = E)
  class(res) <- "SpeciesStats"
  res
}
