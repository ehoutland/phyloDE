#' @title Empirical Bayes estimation of gene-wise phylogenetic parameters
#'
#' @description
#' For every gene, the REML (default) or ML profile log-likelihood of the
#' phylogenetic linear model is evaluated on a fixed grid of correlation
#' structures (Pagel's \eqn{\lambda} for the BM with measurement error, or
#' \eqn{(\rho', \lambda)} for the OU with measurement error, see Details).
#' A non-parametric prior on the grid is then estimated across all genes by
#' maximising the marginal likelihood (with \code{\link[mixsqp]{mixsqp}}),
#' and each gene gets its posterior distribution over the grid.
#' Posterior summaries give regularized, gene-specific phylogenetic parameters:
#' a middle ground between the consensus tree of
#' \code{\link{phylogeneticCorrelations}} (complete pooling) and the
#' independent gene-wise fits of \code{\link{phylolmFit}} with
#' \code{use_consensus = FALSE} (no pooling).
#'
#' @inheritParams phylolmFit
#' @param n_lambda number of grid points for \eqn{\lambda} (ignored if
#' \code{lambda_grid} is given). Default to 41.
#' @param n_rho number of grid points for \eqn{\rho'} (OU only, ignored if
#' \code{rho_grid} is given). Default to 15.
#' @param lambda_grid optional vector of grid values for \eqn{\lambda}, in \eqn{[0, 1)}.
#' @param rho_grid optional vector of grid values for \eqn{\rho' = 1 - \rho}, in \eqn{(0, 1]}.
#' @param keep_posterior if \code{TRUE} (the default), the genes by grid matrix
#' of posterior probabilities is kept in the result. Set to \code{FALSE} to
#' save memory on large data sets.
#' @param prior_smoothing bandwidths (standard deviations, on the logit scale)
#' of the Gaussian bumps used as prior components (see Details).
#' Default to \code{c(0.5, 1, 2)}. Set to \code{NULL} for a non-smoothed prior
#' (free weights on the grid points, the NPMLE).
#' @param prior_pseudocount a small mass added uniformly to the estimated prior
#' so that no grid point has exactly zero prior probability. Default to 1e-8.
#'
#' @details
#' The model for each gene is \eqn{y = X\beta + \epsilon}, \eqn{\epsilon \sim N(0, \sigma^2 C(\theta))},
#' where \eqn{C(\theta)} is the correlation matrix induced by the (replicate)
#' tree transformed by the phylogenetic parameters \eqn{\theta} and re-scaled
#' to unit height. \eqn{\beta} and \eqn{\sigma^2} are profiled out analytically,
#' so that evaluating the likelihood of all genes at one grid point costs a
#' single Cholesky factorisation and a whitening of the data matrix.
#'
#' \itemize{
#' \item For \code{model = "BM"} with \code{measurement_error = TRUE}
#' (or \code{model = "lambda"} without measurement error), \eqn{\theta = \lambda},
#' the proportion of the variance that is phylogenetic, \eqn{C = \lambda V + (1 - \lambda) I}.
#' \item For \code{model = "OUfixedRoot"} with \code{measurement_error = TRUE},
#' \eqn{\theta = (\rho', \lambda)}: the tree is first transformed according to an
#' OU with selection strength \eqn{\alpha}, parametrized by
#' \eqn{\rho' = 1 - \rho = (1 - e^{-2 \alpha h}) / (2 \alpha h)}
#' (the fraction of the variance that can be explained by neutral BM evolution,
#' see \code{\link{rhoFromAlpha}}), and \eqn{\lambda} is then applied to the
#' OU-transformed tree, as in \code{\link{phylolmFit}}.
#' }
#'
#' The prior is a mixture distribution on the grid, estimated by maximum
#' marginal likelihood, exactly as in the \code{ashr} package but using the
#' exact profile likelihood of the phylogenetic model rather than a normal
#' approximation of a point estimate (which is not available for these
#' bounded parameters). As in \code{ashr}, the mixture components overlap so
#' that the estimated prior is smooth: each component is a Gaussian bump on
#' the logit scale of the parameters, centered on a grid point, with a bandwidth
#' taken in \code{prior_smoothing}; the bumps are discretised on the grid.
#' Without smoothing (\code{prior_smoothing = NULL}), the prior is the
#' non-parametric maximum likelihood estimate, which is typically supported on a
#' handful of grid points and gives over-confident gene-wise posteriors.
#' The per-gene posterior mean of \eqn{(\lambda, \rho')}
#' is "snapped" to the nearest grid point, so that genes sharing a grid point
#' share the same correlation structure, which makes the downstream fit
#' (see \code{\link{phylolmFit}}, argument \code{eb_model}) fast.
#'
#' @return An object of class \code{\link{PhyloEBModel-class}}, a list with components
#' \describe{
#' \item{grid}{a data frame with the grid of parameters (\code{lambda}, \code{rho_prime}, \code{rho}, \code{alpha}).}
#' \item{prior}{the estimated prior probabilities on the grid.}
#' \item{posterior}{(if \code{keep_posterior}) the genes by grid matrix of posterior probabilities.}
#' \item{lambda, rho_prime, rho, alpha}{the gene-wise posterior means.}
#' \item{cell}{the index of the grid point nearest to the gene-wise posterior mean.}
#' \item{map}{the index of the gene-wise maximum a posteriori grid point.}
#' \item{logml}{the gene-wise log marginal likelihood under the fitted prior.}
#' \item{loglik_max}{the gene-wise maximum of the profile log-likelihood over the grid.}
#' }
#'
#' @examples
#' data(crayfish)
#' norm_data <- lengthNormalizeRNASeq(crayfish$counts, crayfish$lengths)
#' design <- model.matrix(~ sights, model.frame(crayfish$sights))
#' eb <- phyloEB(norm_data[1:200, ], design = design, phy = crayfish$tree, model = "BM")
#' eb
#' head(getParameters(eb, consensus = FALSE))
#' pfit <- phylolmFit(norm_data[1:200, ], design = design, phy = crayfish$tree,
#'                    model = "BM", eb_model = eb)
#' pfit <- limma::eBayes(pfit)
#' limma::topTable(pfit, coef = 2)
#'
#' @seealso \code{\link{phylolmFit}}, \code{\link{phylogeneticCorrelations}}, \code{\link[mixsqp]{mixsqp}}
#'
#' @export
#'
phyloEB <- function(object, design = NULL, phy, col_species = NULL,
                    model = c("BM", "OUfixedRoot", "lambda"),
                    measurement_error = TRUE,
                    REML = TRUE,
                    n_lambda = 41, n_rho = 15,
                    lambda_grid = NULL, rho_grid = NULL,
                    keep_posterior = TRUE,
                    prior_smoothing = c(0.5, 1, 2),
                    prior_pseudocount = 1e-8) {

  ## Checks (same as phylolmFit)
  y <- check_expression_matrix(object)
  phy <- check_tree(phy, y, col_species)
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
  design <- check_design_matrix(design, y, phy)
  model <- match.arg(model)
  if (model == "BM" && !measurement_error) stop("A BM without measurement error has no phylogenetic parameter to estimate.")
  if (model == "lambda" && measurement_error) stop("Measurement error is not allowed with lambda model.")

  ## Grid
  grid <- phyloEBGrid(phy, model, measurement_error, n_lambda, n_rho, lambda_grid, rho_grid)

  ## Gene by grid profile log-likelihood
  LL <- phyloEBLogLik(y_data, design, phy, model, grid, REML = REML)

  ## Prior: maximum marginal likelihood over mixtures of (smoothed) grid components
  ll_max <- apply(LL, 1, max)
  lik <- exp(LL - ll_max)
  B <- eb_prior_components(grid, prior_smoothing)
  fit <- mixsqp::mixsqp(lik %*% B, control = list(verbose = FALSE))
  weights <- pmax(fit$x, 0)
  weights <- weights / sum(weights)
  prior <- as.vector(B %*% weights)
  prior <- (prior + prior_pseudocount) / (1 + prior_pseudocount * length(prior))

  ## Posterior
  post <- sweep(lik, 2, prior, "*")
  logml <- ll_max + log(rowSums(post))
  post <- post / rowSums(post)

  lambda_mean <- as.vector(post %*% grid$lambda)
  rho_prime_mean <- as.vector(post %*% grid$rho_prime)
  cell <- snap_to_grid(lambda_mean, rho_prime_mean, grid)
  map <- apply(post, 1, which.max)

  res <- list(grid = grid,
              prior = prior,
              posterior = if (keep_posterior) post else NULL,
              lambda = lambda_mean,
              rho_prime = rho_prime_mean,
              rho = 1 - rho_prime_mean,
              alpha = grid$alpha[cell],
              cell = cell,
              map = map,
              logml = logml,
              loglik_max = ll_max,
              model = model,
              measurement_error = measurement_error,
              REML = REML,
              phy = phy,
              design = design,
              genes = rownames(y_data))
  return(new("PhyloEBModel", res))
}

#' @title Prior mixture components on the grid
#'
#' @description
#' Matrix of the discretised prior components: one column per (grid point,
#' bandwidth) pair, giving the probability of each grid point under a Gaussian
#' bump (on the logit scale of \eqn{\lambda} and \eqn{\rho'}) centered on that
#' grid point. Without smoothing, the identity matrix (point masses).
#'
#' @param grid a grid as returned by \code{\link{phyloEBGrid}}.
#' @param bandwidths vector of bandwidths (logit scale), or \code{NULL}.
#'
#' @return A grid points by components matrix, with columns summing to one.
#'
#' @keywords internal
#'
eb_prior_components <- function(grid, bandwidths = NULL) {
  K <- nrow(grid)
  if (is.null(bandwidths) || all(bandwidths <= 0)) return(diag(K))
  logit_clamp <- function(x) stats::qlogis(pmin(pmax(x, 1e-6), 1 - 1e-6))
  d2 <- outer(logit_clamp(grid$lambda), logit_clamp(grid$lambda), "-")^2
  if (length(unique(grid$rho_prime)) > 1) d2 <- d2 + outer(logit_clamp(grid$rho_prime), logit_clamp(grid$rho_prime), "-")^2
  B <- do.call(cbind, lapply(bandwidths, function(s) {
    b <- exp(-0.5 * d2 / s^2)
    sweep(b, 2, colSums(b), "/")
  }))
  return(B)
}

#' @title Grid of phylogenetic parameters
#'
#' @description
#' Build the grid of phylogenetic parameters used by \code{\link{phyloEB}}.
#'
#' @inheritParams phyloEB
#'
#' @return A data frame with one row per grid point and columns
#' \code{lambda}, \code{rho_prime}, \code{rho} and \code{alpha}
#' (\code{alpha = 0}, \code{rho = 0} for the BM).
#'
#' @keywords internal
#'
phyloEBGrid <- function(phy, model, measurement_error, n_lambda = 41, n_rho = 15,
                        lambda_grid = NULL, rho_grid = NULL) {
  if (is.null(lambda_grid)) {
    if (measurement_error || model == "lambda") {
      lambda_max <- getMaxLambda(getMinError(phy))
      lambda_grid <- c(0, stats::plogis(seq(-7, 7, length.out = n_lambda - 2)), lambda_max)
    } else {
      lambda_grid <- 1
    }
  }
  if (any(lambda_grid < 0 | lambda_grid > 1)) stop("`lambda_grid` must be in [0, 1].")
  if (model == "OUfixedRoot") {
    h_tree <- tree_height(phy)
    alpha_bounds <- getBoundsSelectionStrength(phy)
    if (is.null(rho_grid)) {
      rho_prime_bounds <- rho_prime(rev(alpha_bounds), h_tree)
      rho_grid <- stats::plogis(seq(stats::qlogis(rho_prime_bounds[1]) + 0.05, stats::qlogis(rho_prime_bounds[2]) - 0.05, length.out = n_rho))
    }
    if (any(rho_grid <= 0 | rho_grid > 1)) stop("`rho_grid` must be in (0, 1].")
    alpha_grid <- sapply(rho_grid, function(rr) rho_prime_inv(rr, h_tree, alpha_bounds))
    grid <- expand.grid(lambda = lambda_grid, rho_prime = rho_grid)
    grid$alpha <- alpha_grid[match(grid$rho_prime, rho_grid)]
  } else {
    grid <- data.frame(lambda = lambda_grid, rho_prime = 1, alpha = 0)
  }
  grid$rho <- 1 - grid$rho_prime
  rownames(grid) <- NULL
  return(grid)
}

#' @title Transformed tree for one grid point
#'
#' @description
#' Transform the tree according to the phylogenetic parameters of one
#' grid point (OU transform, then lambda transform), and rescale it to unit height.
#'
#' @param phy the (replicate) tree.
#' @param model the phylogenetic model.
#' @param lambda the lambda parameter.
#' @param alpha the OU selection strength (ignored for the BM and lambda models).
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
eb_transform_tree <- function(phy, model, lambda, alpha = 0) {
  tree <- phy
  if (model == "OUfixedRoot") {
    tree <- phylolm::transf.branch.lengths(tree, "OUfixedRoot", parameters = list(alpha = alpha))$tree
  }
  if (lambda < 1) {
    tree <- phylolm::transf.branch.lengths(tree, "lambda", parameters = list(lambda = lambda))$tree
  }
  tree <- rescale_tree(tree)
  return(tree)
}

#' @title Cholesky factors of the grid correlation matrices
#'
#' @description
#' Upper Cholesky factor \eqn{R} (\eqn{C = R^T R}) of the correlation matrix
#' of each grid point, and the log-determinant of \eqn{C}.
#'
#' @inheritParams eb_transform_tree
#' @param grid a grid as returned by \code{\link{phyloEBGrid}}.
#' @param cells the rows of the grid for which to compute the factors (default to all).
#'
#' @return A list with components \code{C_chol} (a list of matrices), \code{logdet} (a vector)
#' and \code{trees} (a list of transformed trees).
#'
#' @keywords internal
#'
eb_grid_chol <- function(phy, model, grid, cells = seq_len(nrow(grid))) {
  trees <- lapply(cells, function(k) eb_transform_tree(phy, model, grid$lambda[k], grid$alpha[k]))
  C_chol <- lapply(trees, function(tt) chol(ape::vcv(tt)))
  logdet <- sapply(C_chol, function(R) 2 * sum(log(diag(R))))
  return(list(C_chol = C_chol, logdet = logdet, trees = trees))
}

#' @title Profile log-likelihood of all genes on the grid
#'
#' @description
#' REML (or ML) profile log-likelihood (with \eqn{\beta} and \eqn{\sigma^2}
#' profiled out) of every gene at every grid point.
#'
#' @param y_data genes by samples matrix of normalized expression, columns in the order of the tree tips.
#' @param design samples by coefficients design matrix.
#' @inheritParams eb_grid_chol
#' @inheritParams phyloEB
#'
#' @return A genes by grid points matrix of log-likelihoods.
#'
#' @keywords internal
#'
phyloEBLogLik <- function(y_data, design, phy, model, grid, REML = TRUE) {
  n <- ncol(y_data)
  p <- ncol(design)
  tY <- t(y_data)
  cst <- log(2 * pi) + 1
  LL <- matrix(NA_real_, nrow = nrow(y_data), ncol = nrow(grid))
  for (k in seq_len(nrow(grid))) {
    tree_k <- eb_transform_tree(phy, model, grid$lambda[k], grid$alpha[k])
    R <- chol(ape::vcv(tree_k))
    logdetC <- 2 * sum(log(diag(R)))
    Yt <- backsolve(R, tY, transpose = TRUE)
    Xt <- backsolve(R, design, transpose = TRUE)
    qrX <- qr(Xt)
    RSS <- colSums(qr.resid(qrX, Yt)^2)
    if (REML) {
      logdetXtX <- 2 * sum(log(abs(diag(qr.R(qrX)))))
      LL[, k] <- -0.5 * ((n - p) * (cst + log(RSS / (n - p))) + logdetC + logdetXtX)
    } else {
      LL[, k] <- -0.5 * (n * (cst + log(RSS / n)) + logdetC)
    }
  }
  rownames(LL) <- rownames(y_data)
  return(LL)
}

#' @title Nearest grid point
#'
#' @description
#' Index of the grid point closest to each gene-wise posterior mean
#' (Euclidean distance in the \eqn{(\lambda, \rho')} plane).
#'
#' @param lambda vector of gene-wise lambda values.
#' @param rho_prime vector of gene-wise rho_prime values.
#' @param grid a grid as returned by \code{\link{phyloEBGrid}}.
#'
#' @return An integer vector of grid indices.
#'
#' @keywords internal
#'
snap_to_grid <- function(lambda, rho_prime, grid) {
  lambdas <- sort(unique(grid$lambda))
  rhos <- sort(unique(grid$rho_prime))
  nearest <- function(x, vals) vals[max.col(-abs(outer(x, vals, "-")), ties.method = "first")]
  l_snap <- nearest(lambda, lambdas)
  r_snap <- nearest(rho_prime, rhos)
  match(paste(l_snap, r_snap), paste(grid$lambda, grid$rho_prime))
}

#' @title Fit using limma on EB-regularized gene-wise correlation structures
#'
#' @description
#' Whiten the data and design with the Cholesky factor of the grid cell of each
#' gene and apply \code{\link[limma]{lmFit}}, one call per grid cell.
#'
#' @param y_data genes by samples matrix of normalized expression.
#' @param design samples by coefficients design matrix.
#' @param eb_model a \code{\link{PhyloEBModel-class}} object.
#' @param ... further parameters to be passed to \code{\link[limma]{lmFit}}.
#'
#' @return A list with the fitted components (\code{coefficients}, \code{sigma},
#' \code{stdev.unscaled}, \code{df.residual}, \code{Amean}, \code{qr},
#' \code{cov.coefficients}, \code{pivot}), the gene-wise \code{C_tree} and
#' \code{trans_tree} (lists sharing the matrices of each cell), and the
#' gene-wise \code{design_trans}.
#'
#' @keywords internal
#'
lmFitEB <- function(y_data, design, eb_model, ...) {
  cells <- sort(unique(eb_model$cell))
  chol_cells <- eb_grid_chol(eb_model$phy, eb_model$model, eb_model$grid, cells)
  G <- nrow(y_data)
  p <- ncol(design)
  coefficients <- matrix(NA_real_, G, p, dimnames = list(rownames(y_data), colnames(design)))
  stdev.unscaled <- coefficients
  sigma <- Amean <- rep(NA_real_, G)
  names(sigma) <- names(Amean) <- rownames(y_data)
  qr_list <- cov_list <- pivot_list <- design_trans <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    genes <- which(eb_model$cell == cells[i])
    R <- chol_cells$C_chol[[i]]
    Xt <- transform_design_one_tree(R, design)
    Yt <- t(backsolve(R, t(y_data[genes, , drop = FALSE]), transpose = TRUE))
    dimnames(Yt) <- list(rownames(y_data)[genes], colnames(y_data))
    fit <- limma::lmFit(Yt, design = Xt, ...)
    coefficients[genes, ] <- fit$coefficients
    stdev.unscaled[genes, ] <- fit$stdev.unscaled
    sigma[genes] <- fit$sigma
    Amean[genes] <- fit$Amean
    qr_list[[i]] <- fit$qr
    cov_list[[i]] <- fit$cov.coefficients
    pivot_list[[i]] <- fit$pivot
    design_trans[[i]] <- Xt
  }
  idx <- match(eb_model$cell, cells)
  return(list(coefficients = coefficients,
              sigma = sigma,
              stdev.unscaled = stdev.unscaled,
              df.residual = rep(ncol(y_data) - p, G),
              Amean = Amean,
              qr = qr_list[idx],
              cov.coefficients = cov_list[idx],
              pivot = pivot_list[idx],
              C_tree = chol_cells$C_chol[idx],
              trans_tree = chol_cells$trees[idx],
              design_trans = design_trans[idx]))
}

#' @title Class PhyloEBModel
#'
#' @description
#' S4 class holding the result of \code{\link{phyloEB}}: a list, see the
#' \code{value} section of \code{\link{phyloEB}}.
#'
#' @name PhyloEBModel-class
#' @rdname PhyloEBModel-class
#'
#' @exportClass PhyloEBModel
#'
setClass("PhyloEBModel", representation("list"))

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("getParameters", "PhyloEBModel", function(object, consensus = TRUE) {
  if (consensus) {
    params <- c(lambda = sum(object$prior * object$grid$lambda))
    if (object$model == "OUfixedRoot") params <- c(params, rho = sum(object$prior * object$grid$rho))
    return(params)
  }
  params <- data.frame(lambda = object$lambda, row.names = object$genes)
  if (object$model == "OUfixedRoot") params$rho <- object$rho
  return(params)
})

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("show", "PhyloEBModel", function(object) {
  cat(is(object)[[1]], "\n",
      "  Empirical Bayes grid on: ", length(object$genes), " genes.\n",
      "  Model: ", object$model, ifelse(object$measurement_error, ", with", ", without"), " measurement error.\n",
      "  Grid: ", length(unique(object$grid$lambda)), " lambda values",
      ifelse(object$model == "OUfixedRoot", paste0(" x ", length(unique(object$grid$rho_prime)), " rho values"), ""), ".\n",
      "  Prior mean: ",
      sep = ""
  )
  getParametersString(object)
  cat("\n")
})

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("plotParameters", "PhyloEBModel", function(object, ...) {
  plot_params_internal(getParameters(object, consensus = TRUE),
                       getParameters(object, consensus = FALSE), ...)
})

check.eb_model <- function(eb_model, model, measurement_error, y_data) {
  if (!inherits(eb_model, "PhyloEBModel")) stop("`eb_model` must be of class `PhyloEBModel`, see function `phyloEB`.")
  if (eb_model$model != model) {
    stop(paste0("The EB model was computed with the ", eb_model$model, " model, but the model used is ", model, ". Please check which model is the correct one."))
  }
  if (eb_model$measurement_error != measurement_error) {
    stop(paste0("The EB model was computed with `measurement_error=", eb_model$measurement_error, "`, but you specified `measurement_error=", measurement_error, "`. Please check which model is the correct one."))
  }
  if (length(eb_model$genes) != nrow(y_data) || !all(eb_model$genes == rownames(y_data))) {
    stop("The genes of the EB model do not match the rows of the expression matrix (same genes, same order).")
  }
}
