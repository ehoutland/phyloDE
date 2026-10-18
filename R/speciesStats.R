#' @title Species-level summary statistics of expression
#'
#' @description
#' Reduce a samples-level data set to one estimate of (log2) abundance per
#' gene and species, with its standard error, the "sufficient statistics"
#' used by \code{\link{phyloMash}} and \code{\link{phyloShiftGLS}}.
#' Two routes are available:
#' \itemize{
#' \item \code{method = "limma"}: the species means of the normalized
#' (log2) data, with a pooled within-species variance per gene, moderated across
#' genes with \code{\link[limma]{squeezeVar}} (the original \code{phyloMash} route).
#' \item \code{method = "edgeR"}: a negative binomial quasi-likelihood fit of
#' the counts with \code{\link[edgeR]{glmQLFit}} and the design
#' \code{~ 0 + species} (plus optional sample-level covariates). The species
#' estimates are the fitted coefficients (converted to log2), and their standard
#' errors come from the negative binomial Fisher information at the fitted means,
#' scaled by the squeezed quasi-likelihood dispersion. The uncertainty of a
#' species thus depends on its own read depth, and the mean-variance
#' relationship is the negative binomial one.
#' }
#'
#' @param object for \code{"limma"}, a matrix of normalized log2 expression
#' (genes by samples, see \code{\link{lengthNormalizeRNASeq}}); for \code{"edgeR"},
#' a matrix of counts or a \code{\link[edgeR]{DGEList}} (with its normalization
#' factors and, if relevant, offsets).
#' @param col_species a character vector giving the species of each sample
#' (column). If \code{NULL}, parsed from the column names (see \code{\link{phylolmFit}}).
#' @param method \code{"limma"} (the default) or \code{"edgeR"}, see Description.
#' @param covariates an optional samples by covariates matrix of sample-level
#' covariates (no intercept) to be regressed out. With \code{"limma"}, the
#' covariates are projected out of the data before taking the species means.
#' @param trend for \code{"limma"}, whether the variance moderation allows for a
#' mean-variance trend (default \code{TRUE}).
#' @param legacy for \code{"edgeR"}, passed to \code{\link[edgeR]{glmQLFit}};
#' \code{FALSE} (the default) uses the new quasi-likelihood method of edgeR 4
#' with adjusted degrees of freedom.
#' @param prior_count for \code{"edgeR"}, the prior count of \code{\link[edgeR]{glmQLFit}}.
#' Species without any read in a gene would otherwise get an estimate of
#' \code{-Inf} with an infinite standard error; with a prior count their
#' estimate and standard error (evaluated at the shrunk means) are finite
#' but are set by the prior count rather than by the data (see \code{dropout}).
#' @param sample_weights an optional vector of sample quality weights (one per
#' sample, e.g. from \code{\link{sampleQualityWeights}}): the variance of a
#' sample's observations is taken as proportional to \code{1 / weight}. With
#' \code{"limma"} they are the array weights of \code{\link[limma]{lmFit}}; with
#' \code{"edgeR"} they are observation weights of \code{\link[edgeR]{glmQLFit}},
#' and the Fisher information is weighted accordingly.
#' @param ... further arguments passed to \code{\link[edgeR]{glmQLFit}}.
#'
#' @return An object of class \code{SpeciesStats}, a list with
#' \describe{
#' \item{mean}{genes by species matrix of log2 abundance estimates.}
#' \item{se}{genes by species matrix of their standard errors.}
#' \item{df}{the degrees of freedom of the t distribution of the estimates
#' (residual plus prior degrees of freedom, a vector for \code{"edgeR"}).}
#' \item{dropout}{genes by species logical matrix: \code{TRUE} if the species
#' has no read at all for the gene (\code{"edgeR"}) or if its mean is the
#' minimum possible (all samples at the pseudo-count floor; \code{"limma"}).}
#' \item{n_rep}{the number of samples of each species.}
#' \item{method, species, genes}{bookkeeping.}
#' }
#'
#' @seealso \code{\link{phyloMash}}, \code{\link{phyloShiftGLS}}
#'
#' @export
#'
speciesStats <- function(object, col_species = NULL, method = c("limma", "edgeR"),
                         covariates = NULL, sample_weights = NULL, trend = TRUE, legacy = FALSE, prior_count = 0.125, ...) {
  method <- match.arg(method)
  if (method == "edgeR") {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is needed for `method = 'edgeR'`.")
    y <- if (inherits(object, "DGEList")) object else edgeR::normLibSizes(edgeR::DGEList(as.matrix(object)))
    counts <- y$counts
  } else {
    counts <- NULL
    y <- check_expression_matrix(object)$exprs
  }
  samples <- if (method == "edgeR") colnames(counts) else colnames(y)
  if (is.null(col_species)) col_species <- sub("(_|\\.).*$", "", samples)
  sp <- factor(col_species)
  S <- nlevels(sp); r <- as.vector(table(sp)); n <- length(sp)
  if (any(r == 0)) stop("Empty species level.")
  X <- stats::model.matrix(~ 0 + sp); colnames(X) <- levels(sp)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != n) stop("`covariates` must have one row per sample.")
    ## covariates orthogonalized to the species space: the species coefficients stay species means
    X <- cbind(X, qr.resid(qr(X), covariates))
  }
  if (!is.null(sample_weights) && (length(sample_weights) != n || any(sample_weights <= 0))) stop("`sample_weights` must be one positive weight per sample.")

  if (method == "limma") {
    G <- nrow(y)
    df_res <- n - ncol(X)
    if (df_res < 2) stop("Too few residual degrees of freedom to estimate the within-species variance.")
    fit <- limma::lmFit(y, X, weights = sample_weights)
    M <- fit$coefficients[, seq_len(S), drop = FALSE]
    sq <- if (trend) limma::squeezeVar(fit$sigma^2, df = df_res, covariate = fit$Amean) else limma::squeezeVar(fit$sigma^2, df = df_res)
    s2post <- sq$var.post
    df <- df_res + if (is.finite(sq$df.prior)) sq$df.prior else 0
    se <- sqrt(s2post) * fit$stdev.unscaled[, seq_len(S), drop = FALSE]
    dropout <- M <= rep(apply(y, 1, min), S) + 1e-8
    dimnames(se) <- dimnames(dropout) <- dimnames(M)
    if (any(r == 1)) message("Species with a single sample: ", paste(levels(sp)[r == 1], collapse = ", "), ". Their within-species variance is pooled from the other species.")
  } else {
    if (!is.null(sample_weights)) y$weights <- matrix(sample_weights, nrow(counts), n, byrow = TRUE)
    if (is.null(y$common.dispersion) && is.null(y$trended.dispersion) && is.null(y$tagwise.dispersion) && legacy) y <- edgeR::estimateDisp(y, X)
    fit <- edgeR::glmQLFit(y, X, legacy = legacy, prior.count = prior_count, ...)
    G <- nrow(counts)
    beta <- fit$coefficients[, seq_len(S), drop = FALSE]
    ## NB Fisher information per species (diagonal: species are disjoint blocks, up to the covariates,
    ## whose effect on the species blocks is ignored), evaluated at the means implied by the reported
    ## (prior-count shrunk) coefficients: for a species without any read, the fitted mean is 0 and the
    ## information would vanish, while the coefficient itself is finite only because of the prior count.
    mu <- exp(fit$coefficients %*% t(X) + as.matrix(fit$offset))
    phi <- fit$dispersion
    if (length(phi) == 1) phi <- rep(phi, G)
    wq <- if (is.null(sample_weights)) 1 else rep(sample_weights, each = G)
    info <- (wq * mu / (1 + phi * mu)) %*% X[, seq_len(S), drop = FALSE]
    var_post <- if (!is.null(fit$var.post)) fit$var.post else fit$s2.post
    se <- sqrt(var_post / info)
    M <- beta / log(2); se <- se / log(2)
    df_res <- if (!is.null(fit$df.residual.adj)) fit$df.residual.adj else fit$df.residual
    df <- df_res + fit$df.prior
    dropout <- (counts %*% X[, seq_len(S), drop = FALSE]) == 0
    dimnames(M) <- dimnames(se) <- dimnames(dropout) <- list(rownames(counts), levels(sp))
  }
  res <- list(mean = M, se = se, df = df, dropout = dropout, n_rep = setNames(r, levels(sp)),
              method = method, species = levels(sp), genes = rownames(M))
  class(res) <- "SpeciesStats"
  res
}

#' @title Phylogenetic generalized least squares on species-level statistics
#'
#' @description
#' Fit, for every gene, the phylogenetic linear model at the species level,
#' \eqn{\hat{m}_g \sim N(X \beta_g, \sigma^2_g C(\theta_g) + D_g)}, where
#' \eqn{\hat{m}_g} are the species-level estimates of \code{\link{speciesStats}},
#' \eqn{D_g} the diagonal matrix of their squared standard errors (known),
#' \eqn{C(\theta)} the phylogenetic correlation among species and
#' \eqn{\sigma^2_g} the phylogenetic (between-species) variance. The
#' phylogenetic parameters \eqn{\theta} (\eqn{\rho'} for the OU; the BM has
#' none) and the ratio of phylogenetic to measurement variance are handled on
#' a grid with the empirical Bayes machinery of \code{\link{phyloEB}}
#' (profile likelihood on the grid, smooth prior across genes, posterior
#' averaging of the estimates), so that the result is a species-level analogue
#' of \code{\link{phyloShiftFit}} with heteroscedastic, read-depth-aware
#' measurement error: with \code{method = "edgeR"} statistics, a species with
#' few reads weighs less than one with many.
#'
#' Two kinds of tests are available, through \code{design} and \code{shifts}:
#' a coefficient of a species-level design (e.g. a condition such as
#' blind/sighted), and/or branch-specific shifts as in \code{\link{phyloShiftFit}}.
#'
#' @param stats a \code{SpeciesStats} object.
#' @param phy the species tree (one tip per species), ultrametric.
#' @param design an optional species by coefficients design matrix (with an
#' intercept); the last coefficient is tested. If \code{NULL}, only the intercept.
#' @param shifts if \code{TRUE} (the default), also test a shift on every
#' testable branch of the tree (indicator added to \code{design}, one at a time).
#' @param model \code{"BM"} (the default) or \code{"OUfixedRoot"}.
#' @param n_lambda,n_rho number of grid values for the phylogenetic variance
#' (relative to the gene's mean squared standard error, log-spaced from 0.05 to 50,
#' plus zero) and for \eqn{\rho'}. The reported \code{lambda} is the share of
#' phylogenetic variance, \eqn{\tau^2 / (1 + \tau^2)}.
#' @param prior_smoothing,prior_pseudocount see \code{\link{phyloEB}}.
#'
#' @details
#' For a given grid point, the covariance is \eqn{\bar{s}^2 (\tau^2 C + D / \bar{s}^2)},
#' with \eqn{\bar{s}^2} the mean squared standard error of the gene: the scale
#' is known (from the measurement error), not profiled, which is what makes the
#' standard errors of the coefficients valid with as few as four species.
#' The REML log-likelihood is evaluated in closed form for every gene and grid
#' point (one Cholesky per gene and grid point, cheap because the matrices are
#' species by species), the prior over the grid is estimated across genes as
#' in \code{\link{phyloEB}}, and the estimates are averaged over each gene's posterior.
#'
#' @return An object of class \code{PhyloShiftFit} (see \code{\link{phyloShiftFit}}),
#' with the design coefficient as an extra column named after it, if \code{design}
#' is given, in \code{Bhat}, \code{Shat}, \code{t} and \code{p.value}.
#'
#' @seealso \code{\link{speciesStats}}, \code{\link{phyloShiftFit}}, \code{\link{phyloMash}}
#'
#' @export
#'
phyloShiftGLS <- function(stats, phy, design = NULL, shifts = TRUE, model = c("BM", "OUfixedRoot"),
                          n_lambda = 21, n_rho = 9, prior_smoothing = c(0.5, 1, 2), prior_pseudocount = 1e-8) {
  if (!inherits(stats, "SpeciesStats")) stop("`stats` must be a `SpeciesStats` object.")
  model <- match.arg(model)
  if (!all(stats$species %in% phy$tip.label)) stop("Some species of `stats` are not in the tree.")
  phy <- ape::keep.tip(phy, stats$species)
  ord <- match(phy$tip.label, stats$species)
  M <- stats$mean[, ord, drop = FALSE]; SE <- stats$se[, ord, drop = FALSE]
  S <- ncol(M); G <- nrow(M)
  df_meas <- stats$df
  if (is.null(design)) design <- matrix(1, S, 1, dimnames = list(phy$tip.label, "(Intercept)"))
  design <- as.matrix(design)
  if (nrow(design) != S) stop("`design` must have one row per species.")
  if (!is.null(rownames(design))) design <- design[phy$tip.label, , drop = FALSE]
  tests <- list()
  if (ncol(design) > 1) tests[[colnames(design)[ncol(design)]]] <- list(X0 = design[, -ncol(design), drop = FALSE], z = design[, ncol(design)])
  if (shifts) {
    bi <- branchIndicators(phy, phy$tip.label)
    keep <- bi$branches$testable
    ## drop branches whose indicator is collinear with the design (e.g. a condition that is a clade)
    estimable <- sapply(seq_len(nrow(bi$branches)), function(b) qr(cbind(design, bi$Z[, b]))$rank > ncol(design))
    if (any(keep & !estimable)) message("Branches not estimable given the design, skipped: ", paste(bi$branches$label[keep & !estimable], collapse = ", "))
    keep <- keep & estimable
    for (b in which(keep)) tests[[bi$branches$label[b]]] <- list(X0 = design, z = bi$Z[, b])
    branches <- bi$branches[keep, ]
    Z <- bi$Z[, keep, drop = FALSE]
  } else { branches <- NULL; Z <- NULL }

  ## grid: tau2 = phylogenetic variance in units of the gene's mean squared SE (log scale); rho' for the OU
  tr_u <- rescale_tree(phy)
  tau_grid <- c(0, exp(seq(log(0.05), log(50), length.out = n_lambda - 1)))
  if (model == "OUfixedRoot") {
    ab <- getBoundsSelectionStrength(tr_u)
    rb <- rho_prime(rev(ab), 1); rb <- c(max(rb[1], 1e-3), min(rb[2], 1 - 1e-3))
    rho_grid <- stats::plogis(seq(stats::qlogis(rb[1]) + 0.05, stats::qlogis(rb[2]) - 0.05, length.out = n_rho))
    grid <- expand.grid(tau2 = tau_grid, rho_prime = rho_grid)
    grid$alpha <- sapply(grid$rho_prime, function(r) rho_prime_inv(r, 1, ab))
  } else grid <- data.frame(tau2 = tau_grid, rho_prime = 1, alpha = 0)
  grid$rho <- 1 - grid$rho_prime
  grid$lambda <- grid$tau2 / (1 + grid$tau2)
  K <- nrow(grid)
  Cs <- lapply(seq_len(K), function(k) ape::vcv(eb_transform_tree(tr_u, model, 1, grid$alpha[k])))
  s2bar <- rowMeans(SE^2)
  Dn <- SE^2 / s2bar   # normalized measurement variances (mean 1 per gene)

  ## REML log-likelihood with known scale (no profiling of sigma^2) and FWL estimate, gene by gene
  run_test <- function(X0, z) {
    p <- ncol(X0) + 1
    LL <- est <- se2 <- matrix(NA_real_, G, K)
    for (k in seq_len(K)) {
      Ck <- Cs[[k]]; tau2 <- grid$tau2[k]
      for (g in seq_len(G)) {
        Sig <- (tau2 * Ck + diag(Dn[g, ], S)) * s2bar[g]
        R <- chol(Sig)
        X0t <- backsolve(R, X0, transpose = TRUE); zt <- backsolve(R, z, transpose = TRUE); yt <- backsolve(R, M[g, ], transpose = TRUE)
        q0 <- qr(X0t); zp <- qr.resid(q0, zt); yp <- qr.resid(q0, yt)
        nzk <- sum(zp^2); b <- sum(yp * zp) / nzk
        RSS <- sum(yp^2) - b^2 * nzk
        est[g, k] <- b; se2[g, k] <- 1 / nzk
        LL[g, k] <- -0.5 * (RSS + (S - p) * log(2 * pi) + 2 * sum(log(diag(R))) + 2 * sum(log(abs(diag(qr.R(qr(cbind(X0t, zt))))))))
      }
    }
    pp <- eb_fit_prior(LL, grid, prior_smoothing, prior_pseudocount)
    W <- pp$posterior
    Bhat <- rowSums(W * est)
    Shat <- sqrt(rowSums(W * (se2 + (est - Bhat)^2)))
    ## the measurement variances are estimated with df_meas degrees of freedom: t reference with those df
    dfm <- if (length(df_meas) == 1) rep(df_meas, G) else df_meas
    pv <- rowSums(W * 2 * stats::pt(-abs(est) / sqrt(se2), df = dfm))
    list(Bhat = Bhat, Shat = Shat, p = pv, lambda = as.vector(W %*% grid$lambda), rho = as.vector(W %*% grid$rho), prior = pp$prior)
  }
  out <- lapply(tests, function(tt) run_test(tt$X0, tt$z))
  nm <- names(tests)
  Bhat <- sapply(out, `[[`, "Bhat"); Shat <- sapply(out, `[[`, "Shat"); pval <- sapply(out, `[[`, "p")
  lambda <- sapply(out, `[[`, "lambda"); rho <- sapply(out, `[[`, "rho")
  if (is.null(dim(Bhat))) { Bhat <- matrix(Bhat, G); Shat <- matrix(Shat, G); pval <- matrix(pval, G); lambda <- matrix(lambda, G); rho <- matrix(rho, G) }
  dimnames(Bhat) <- dimnames(Shat) <- dimnames(pval) <- dimnames(lambda) <- dimnames(rho) <- list(stats$genes, nm)
  res <- list(Bhat = Bhat, Shat = Shat, t = Bhat / Shat, p.value = pval, df = NA,
              lambda = lambda, rho = rho, nuisance = "eb_gls", branches = branches, Z = Z,
              priors = sapply(out, `[[`, "prior"), grid = grid, model = model, phy = phy, stats_method = stats$method,
              genes = stats$genes)
  class(res) <- "PhyloShiftFit"
  res
}

#' @title Sample quality weights
#'
#' @description
#' One precision weight per sample, estimated across genes as in
#' \code{\link[limma]{arrayWeights}} (via \code{\link[limma]{voomWithQualityWeights}}):
#' a sample whose observations scatter more around the fitted values of the
#' species (plus covariates) design than the others gets a weight below 1.
#' Outlier libraries are downweighted rather than removed. The weights can be
#' passed to \code{\link{speciesStats}} (\code{sample_weights}) and
#' \code{\link{phyloWeightedFit}} (\code{sample_weights}).
#'
#' @param counts genes by samples matrix of counts, or a \code{\link[edgeR]{DGEList}}.
#' @param col_species the species of each sample.
#' @param covariates optional samples by covariates matrix (no intercept), as in \code{\link{speciesStats}}.
#' @param norm.factors optional normalization factors (TMM by default).
#' @param ... further arguments passed to \code{\link[limma]{voomWithQualityWeights}}.
#'
#' @return A named numeric vector of weights (one per sample).
#'
#' @export
#'
sampleQualityWeights <- function(counts, col_species, covariates = NULL, norm.factors = NULL, ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is needed for function 'sampleQualityWeights'.")
  y <- if (inherits(counts, "DGEList")) counts else edgeR::DGEList(as.matrix(counts))
  if (!is.null(norm.factors)) y$samples$norm.factors <- norm.factors else if (all(y$samples$norm.factors == 1)) y <- edgeR::normLibSizes(y)
  sp <- factor(col_species)
  X <- if (nlevels(sp) > 1) stats::model.matrix(~ 0 + sp) else matrix(1, ncol(y), 1)
  colnames(X) <- levels(sp)
  if (!is.null(covariates)) X <- cbind(X, as.matrix(covariates))
  v <- limma::voomWithQualityWeights(y, X, plot = FALSE, ...)
  w <- v$targets$sample.weights
  names(w) <- colnames(y)
  w
}
