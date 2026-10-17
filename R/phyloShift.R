#' @title Branch indicators for shift detection
#'
#' @description
#' Build the samples by branches indicator matrix used to test for an
#' expression shift on each branch of the species tree: \code{Z[i, b] = 1} if
#' sample \code{i} belongs to a species descending from branch \code{b}.
#' All the branches (edges) of the tree are returned, but only the "testable"
#' ones are used in \code{\link{phyloShiftFit}}: when the root has two child
#' edges, a shift on one of them is the same as a shift of opposite sign on the
#' other (given the intercept), so the one leading to the larger clade is
#' marked as not testable.
#'
#' @param phy the species tree (one tip per species).
#' @param col_species a character vector giving the species of each sample.
#'
#' @return A list with \code{Z} (the samples by branches 0/1 matrix) and
#' \code{branches}, a data frame with one row per edge of \code{phy}:
#' \code{edge}, \code{parent}, \code{child} (node numbers), \code{label}
#' (the species name for a terminal branch, \code{nodeX} otherwise),
#' \code{species} (the descendant species, "+" separated), \code{n_species},
#' \code{terminal} and \code{testable}.
#'
#' @export
#'
branchIndicators <- function(phy, col_species) {
  if (!inherits(phy, "phylo")) stop("object 'phy' must be of class 'phylo'.")
  S <- length(phy$tip.label)
  if (!all(col_species %in% phy$tip.label)) stop("Some species in `col_species` are not in the tree.")
  root <- S + 1
  child <- phy$edge[, 2]
  tips_under <- lapply(child, function(v) if (v <= S) v else { d <- phytools::getDescendants(phy, v); d[d <= S] })
  branches <- data.frame(edge = seq_len(nrow(phy$edge)), parent = phy$edge[, 1], child = child,
                         n_species = lengths(tips_under), terminal = child <= S, stringsAsFactors = FALSE)
  branches$species <- sapply(tips_under, function(v) paste(phy$tip.label[v], collapse = "+"))
  branches$label <- ifelse(branches$terminal, phy$tip.label[pmin(child, S)], paste0("node", child))
  root_edges <- which(branches$parent == root)
  branches$testable <- TRUE
  if (length(root_edges) == 2) branches$testable[root_edges[which.max(branches$n_species[root_edges])]] <- FALSE
  Z <- sapply(tips_under, function(v) as.numeric(col_species %in% phy$tip.label[v]))
  Z <- matrix(Z, nrow = length(col_species), dimnames = list(NULL, branches$label))
  return(list(Z = Z, branches = branches))
}

#' @title Branch-wise shift estimates for all genes
#'
#' @description
#' For every gene and every testable branch of the species tree, estimate the
#' shift in expression carried by that branch (an indicator of the descendant
#' samples added to the null design), with the phylogenetic linear model.
#' The gene-wise phylogenetic parameters are treated with the empirical Bayes
#' grid of \code{\link{phyloEB}}: for each branch, the prior is estimated under
#' the alternative of that branch, and the estimate and its standard error are
#' averaged over the gene-wise posterior on the grid (see Details).
#' Estimates, moderated standard errors, the null correlation among the branch
#' estimates and the expected "leakage" of a shift on one branch into the
#' tests of the other branches are returned, ready for
#' \code{\link{phyloShiftMash}} or \code{\link{phyloShiftAsh}}.
#'
#' @inheritParams phylolmFit
#' @param phy the species tree (one tip per species, see \code{\link{branchIndicators}}).
#' @param design the null design matrix (samples by coefficients), typically
#' the intercept (the default) plus sample-level covariates. The branch
#' indicators are added to it one at a time.
#' @param nuisance how the gene-wise phylogenetic parameters are treated:
#' \code{"eb"} (the default) averages each estimate over the gene-wise posterior
#' on the grid, with a prior estimated under the alternative of each branch;
#' \code{"eb_plugin"} uses the single grid point nearest to the posterior mean;
#' \code{"eb_null"} uses the posterior under the null design for all branches
#' (faster but biased for the genes that do carry a shift, see Details);
#' \code{"consensus"} uses, for each branch, one value for all genes, the trimmed
#' mean of the gene-wise maximum likelihood estimates on the grid, as in
#' \code{\link{phylogeneticCorrelations}}.
#' @param squeeze if \code{TRUE} (the default), the gene-wise residual variances
#' are moderated with \code{\link[limma]{squeezeVar}} before computing the
#' standard errors.
#' @param trim trimming proportions for \code{nuisance = "consensus"}, see \code{\link{phylogeneticCorrelations}}.
#' @param ... further arguments passed to \code{\link{phyloEB}} (grid and prior options).
#'
#' @details
#' Each gene is whitened with the Cholesky factor of a correlation structure.
#' With \eqn{\tilde{y}} and \eqn{\tilde{X}_0} the whitened data and null design,
#' and \eqn{\tilde{z}_b} the whitened indicator of branch \eqn{b}, residualized
#' on \eqn{\tilde{X}_0}, the estimate of the shift on branch \eqn{b} is
#' \eqn{\hat{\beta}_b = \tilde{z}_b^T \tilde{y}^\perp / \|\tilde{z}_b\|^2}
#' (Frisch-Waugh-Lovell), which is exactly the last coefficient of the
#' phylogenetic regression on \code{cbind(design, Z[, b])}.
#'
#' The phylogenetic parameters of a gene that carries a shift are badly biased
#' when estimated under the null design (the shift looks like extra
#' phylogenetic variance), which inflates the standard error of the shift
#' and destroys power. For each branch, the grid likelihood and the prior are
#' therefore computed under the alternative of that branch, with
#' \code{cbind(design, Z[, b])} as design. Plugging in a single value of the
#' parameters is in turn anti-conservative for the genes whose posterior is
#' wide (few species give little information on the gene-wise parameters), so
#' with \code{nuisance = "eb"} the estimate is the posterior mean over the grid
#' of the grid-wise estimates, and its variance is the posterior mean of the
#' grid-wise variances plus the posterior variance of the grid-wise estimates
#' (law of total variance); the p-value is the posterior average of the
#' grid-wise t-test p-values.
#'
#' Under the null, the estimates of the different branches of a gene are
#' correlated, with correlation matrix \eqn{V = \mathrm{cor}(\tilde{Z}^\perp)}
#' (one per grid cell of the null fit). A true shift \eqn{\delta} on branch \eqn{b} gives
#' expected estimates \eqn{\delta A[, b]} on all branches, with
#' \eqn{A[b', b] = \tilde{z}_{b'}^T \tilde{z}_b / \|\tilde{z}_{b'}\|^2}: a shift
#' shared by a clade is attenuated in the tests of the clade members (whose
#' phylogenetic prediction borrows from their shifted relatives) and diluted
#' in the tests of the ancestral branches. These patterns are the canonical
#' covariances of \code{\link{phyloShiftMash}}. Both are computed with the
#' correlation structures of the null fit.
#'
#' @return An object of class \code{PhyloShiftFit}, a list with
#' \describe{
#' \item{Bhat, Shat, t, p.value}{genes by testable branches matrices of estimates, (moderated) standard errors, t statistics and two-sided p-values.}
#' \item{df}{the degrees of freedom of the t statistics (residual plus prior degrees of freedom).}
#' \item{lambda, rho}{genes by branches matrices of the posterior mean parameters under each branch's alternative.}
#' \item{cell}{the grid cell of each gene under the null design.}
#' \item{V, A}{lists (one element per null cell) of the null correlation and leakage matrices among branches.}
#' \item{branches}{the branch table (testable branches only).}
#' \item{Z}{the samples by testable branches indicator matrix.}
#' \item{eb_null}{the \code{\link{PhyloEBModel-class}} of the null design.}
#' \item{priors}{the branch-wise priors on the grid.}
#' \item{s2.prior, df.prior}{the variance prior from \code{\link[limma]{squeezeVar}}.}
#' }
#'
#' @seealso \code{\link{phyloEB}}, \code{\link{phyloShiftMash}}, \code{\link{phyloShiftAsh}}
#'
#' @export
#'
phyloShiftFit <- function(object, design = NULL, phy, col_species = NULL,
                          model = c("BM", "OUfixedRoot", "lambda"),
                          measurement_error = TRUE,
                          nuisance = c("eb", "eb_plugin", "eb_null", "consensus"),
                          squeeze = TRUE, trim = c(0.25, 0.05), ...) {
  ## checks
  y <- check_expression_matrix(object)
  if (is.null(col_species)) col_species <- parse_species(phy, colnames(y$exprs))
  phy_rep <- check_tree(phy, y, col_species)
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy_rep)
  design <- check_design_matrix(design, y, phy_rep)
  model <- match.arg(model)
  nuisance <- match.arg(nuisance)
  dots <- list(...)
  prior_smoothing <- if (!is.null(dots$prior_smoothing)) dots$prior_smoothing else c(0.5, 1, 2)
  prior_pseudocount <- if (!is.null(dots$prior_pseudocount)) dots$prior_pseudocount else 1e-8
  REML <- if (!is.null(dots$REML)) dots$REML else TRUE

  ## branch indicators, in the sample order of y_data
  sample_species <- col_species[match(colnames(y_data), colnames(y$exprs))]
  bi <- branchIndicators(phy, sample_species)
  keep <- bi$branches$testable
  Z <- bi$Z[, keep, drop = FALSE]
  branches <- bi$branches[keep, ]
  rownames(Z) <- colnames(y_data)
  R <- ncol(Z)
  n <- ncol(y_data); p0 <- ncol(design)
  df1 <- n - p0 - 1
  G <- nrow(y_data)

  ## null fit: cells for V and A, prior for the variances
  eb_null <- phyloEB(y_data, design = design, phy = phy_rep, model = model,
                     measurement_error = measurement_error, keep_posterior = TRUE, ...)
  grid <- eb_null$grid
  K <- nrow(grid)
  lam_eps <- pmin(grid$lambda, 1 - 1e-12)

  ## variance prior from the null-design residuals (the grid pass of the first branch gives them)
  Bhat <- Shat <- pval <- lambda <- rho <- matrix(NA_real_, G, R, dimnames = list(rownames(y_data), branches$label))
  priors <- matrix(NA_real_, K, R, dimnames = list(NULL, branches$label))
  s2.prior <- NA; df.prior <- 0
  for (b in seq_len(R)) {
    gp <- phyloEBGridPass(y_data, design, Z[, b], phy_rep, model, grid, REML = REML)
    if (b == 1) {
      ## null-design residual variances: RSS of the model with z, plus the part explained by z, at the null posterior mode
      k0 <- eb_null$map
      RSS0 <- (gp$s2[cbind(seq_len(G), k0)] * df1 + gp$est[cbind(seq_len(G), k0)]^2 * gp$nz[k0])
      if (squeeze) {
        sq <- limma::squeezeVar(RSS0 / (n - p0), df = n - p0)
        df.prior <- sq$df.prior; s2.prior <- sq$var.prior
        if (!is.finite(df.prior)) df.prior <- 0
      }
    }
    s2post <- if (squeeze && df.prior > 0) (df.prior * s2.prior + df1 * gp$s2) / (df.prior + df1) else gp$s2
    df <- df1 + df.prior
    se2 <- sweep(s2post, 2, gp$nz, "/")
    ## posterior over the grid under the alternative of branch b
    if (nuisance == "eb_null") {
      W <- eb_null$posterior
      priors[, b] <- eb_null$prior
    } else if (nuisance == "consensus") {
      mle <- apply(gp$loglik, 1, which.max)
      lam_c <- tanh(mean_trim(atanh(lam_eps[mle]), trim = trim, na.rm = TRUE))
      rho_c <- if (model == "OUfixedRoot") tanh(mean_trim(atanh(pmax(-1, grid$rho[mle])), trim = trim, na.rm = TRUE)) else 0
      W <- matrix(0, G, K); W[, snap_to_grid(lam_c, 1 - rho_c, grid)] <- 1
      priors[, b] <- colMeans(W)
    } else {
      pp <- eb_fit_prior(gp$loglik, grid, prior_smoothing, prior_pseudocount)
      W <- pp$posterior
      priors[, b] <- pp$prior
      if (nuisance == "eb_plugin") {
        cell_b <- snap_to_grid(as.vector(W %*% grid$lambda), as.vector(W %*% grid$rho_prime), grid)
        W <- matrix(0, G, K); W[cbind(seq_len(G), cell_b)] <- 1
      }
    }
    lambda[, b] <- as.vector(W %*% grid$lambda); rho[, b] <- as.vector(W %*% grid$rho)
    Bhat[, b] <- rowSums(W * gp$est)
    Shat[, b] <- sqrt(rowSums(W * (se2 + (gp$est - Bhat[, b])^2)))
    pval[, b] <- rowSums(W * 2 * stats::pt(-abs(gp$est) / sqrt(se2), df = df))
  }
  tstat <- Bhat / Shat

  ## null correlation and leakage among branches, from the null cells
  cells0 <- sort(unique(eb_null$cell))
  cell <- match(eb_null$cell, cells0)
  chol0 <- eb_grid_chol(phy_rep, model, grid, cells0)$C_chol
  V <- A <- cell_info <- vector("list", length(cells0))
  for (i in seq_along(cells0)) {
    Rc <- chol0[[i]]
    q0 <- qr(backsolve(Rc, design, transpose = TRUE))
    Zp <- qr.resid(q0, backsolve(Rc, Z, transpose = TRUE))
    ZZ <- crossprod(Zp)
    V[[i]] <- stats::cov2cor(ZZ)
    A[[i]] <- ZZ / colSums(Zp^2)
    dimnames(V[[i]]) <- dimnames(A[[i]]) <- list(branches$label, branches$label)
    cell_info[[i]] <- list(grid_cell = cells0[i], chol = Rc, Zp = Zp, nz = colSums(Zp^2))
  }

  res <- list(Bhat = Bhat, Shat = Shat, t = tstat, p.value = pval, df = df,
              lambda = lambda, rho = rho, nuisance = nuisance,
              cell = cell, V = V, A = A, cell_info = cell_info, branches = branches, Z = Z,
              eb_null = eb_null, priors = priors, model = model, measurement_error = measurement_error,
              s2.prior = s2.prior, df.prior = df.prior, design = design,
              phy = phy, phy_rep = phy_rep, col_species = sample_species, genes = rownames(y_data))
  class(res) <- "PhyloShiftFit"
  return(res)
}

#' @title Empirical Bayes shift detection across branches with mash
#'
#' @description
#' Run \code{\link[mashr]{mash}} on the branch-wise shift estimates of
#' \code{\link{phyloShiftFit}}: the effects are the shifts on all testable
#' branches, and the prior is a mixture of (scaled) covariance patterns
#' learned from all genes. Three kinds of canonical patterns are used:
#' \itemize{
#' \item "shift" patterns, one per branch: a single shift on that branch,
#' including its expected leakage into the other branch tests (see
#' \code{\link{phyloShiftFit}});
#' \item "drift" patterns: the covariance of the branch estimates of a gene
#' without any shift, whose phylogenetic variance is larger than assumed by
#' its standard errors (few species give little information on the gene-wise
#' phylogenetic parameters, so the null covariance of a gene is uncertain;
#' without these patterns, genes with a large phylogenetic variance are absorbed
#' by the shift patterns);
#' \item optionally, the identity (unstructured effects).
#' }
#' A gene carries a shift on branch \eqn{b} with posterior probability equal
#' to the posterior weight of the shift patterns of \eqn{b}
#' (\code{branch_prob}); the null, drift and identity patterns count as "no shift".
#'
#' @param fit a \code{PhyloShiftFit} object.
#' @param pattern_cells which grid cells to use for the shift patterns:
#' \code{"quantiles"} (the default) uses the cells of the genes at the 10, 50 and 90 percent
#' quantiles of the gene-wise lambda, \code{"modal"} the most populated cell only, or a
#' vector of cell indices.
#' @param drift_lambda values of \eqn{\lambda} for the drift patterns (computed
#' with the whitening of the most populated cell). Default to \code{c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95)}.
#' @param drift_rho_prime values of \eqn{\rho'} for the drift patterns (OU model only).
#' @param identity if \code{TRUE}, add the identity pattern. Default to \code{FALSE}.
#' @param singletons if \code{TRUE}, also include the singleton patterns
#' (one branch moves, the others do not). Default to \code{FALSE}.
#' @param data_driven if \code{TRUE}, add data-driven patterns
#' (\code{\link[mashr]{cov_pca}} then \code{\link[mashr]{cov_ed}}) learned
#' from the genes with a significant condition-by-condition signal. They are
#' counted as "no shift" in \code{branch_prob}. Default to \code{FALSE}.
#' @param n_pca number of principal components for the data-driven patterns.
#' @param lfsr_strong lfsr threshold defining the strong genes for the data-driven patterns.
#' @param null_cor how the correlation among the branch estimates of a gene
#' under the null is obtained: \code{"empirical"} (the default) estimates it from
#' the genes with weak signal on all branches with
#' \code{\link[mashr]{estimate_null_correlation_simple}}; \code{"model"} uses the
#' model-based correlation of the null fit (\code{fit$V}), averaged over the
#' grid cells weighted by their number of genes (\code{shared_V = TRUE}), or gene
#' by gene (\code{shared_V = FALSE}, which forces the slow pure R implementation of \code{mash}).
#' @param shared_V see \code{null_cor}.
#' @param ridge with the same correlation structure for all branches, the branch
#' estimates of a gene are linearly dependent (the estimate of an internal branch is
#' a combination of those of its child branches), so their null correlation matrix is
#' singular or nearly so. It is replaced by \code{(1 - ridge) * V + ridge * I}
#' (default \code{ridge = 0.05}) for \code{mash}.
#' @param ... further arguments passed to \code{\link[mashr]{mash}}.
#'
#' @return A list with the \code{mash} fit (\code{mash}), the genes by branches
#' matrices \code{lfsr} (local false sign rate of the branch effects, which counts
#' drift as an effect), \code{pm} (posterior means) and \code{psd} (posterior standard deviations),
#' \code{branch_prob} (posterior probability that the gene carries a single shift on each branch),
#' \code{shift_prob} (posterior probability of a shift on any branch),
#' \code{top_branch} (the branch with the highest posterior probability),
#' \code{pattern_weights} (the estimated prior weight of each pattern, summed over scales)
#' and \code{V} (the null correlation used).
#'
#' @seealso \code{\link{phyloShiftFit}}, \code{\link{phyloShiftAsh}}, \code{\link[mashr]{mash}}
#'
#' @export
#'
phyloShiftMash <- function(fit, pattern_cells = "quantiles",
                           drift_lambda = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95),
                           drift_rho_prime = c(0.1, 0.5, 0.9),
                           identity = FALSE, singletons = FALSE,
                           data_driven = FALSE, n_pca = 3, lfsr_strong = 0.05,
                           null_cor = c("empirical", "model"), shared_V = TRUE,
                           ridge = 0.05, ...) {
  if (!requireNamespace("mashr", quietly = TRUE)) stop("Package 'mashr' is needed for function 'phyloShiftMash'.")
  R <- ncol(fit$Bhat); G <- nrow(fit$Bhat)
  labels <- fit$branches$label
  modal <- as.integer(names(which.max(table(fit$cell))))
  ## null correlation among the branch estimates
  null_cor <- match.arg(null_cor)
  dots <- list(...)
  regularize <- function(V) (1 - ridge) * V + ridge * diag(R)
  if (null_cor == "empirical") {
    data0 <- mashr::mash_set_data(fit$Bhat, fit$Shat, df = fit$df)
    Vdat <- regularize(mashr::estimate_null_correlation_simple(data0))
  } else if (shared_V || length(fit$V) == 1) {
    wcell <- tabulate(fit$cell, nbins = length(fit$V)) / G
    Vdat <- regularize(Reduce(`+`, Map(`*`, fit$V, wcell)))
  } else {
    Vdat <- array(unlist(lapply(fit$V[fit$cell], regularize)), dim = c(R, R, G))
    if (is.null(dots$algorithm.version)) dots$algorithm.version <- "R"
  }
  data <- mashr::mash_set_data(fit$Bhat, fit$Shat, df = fit$df, V = Vdat)

  ## shift patterns: leakage of a single shift
  if (identical(pattern_cells, "quantiles")) {
    lam <- fit$eb_null$lambda
    cells_p <- unique(fit$cell[order(lam)][pmax(1, round(c(0.1, 0.5, 0.9) * G))])
  } else if (identical(pattern_cells, "modal")) {
    cells_p <- modal
  } else cells_p <- pattern_cells
  Ulist <- list()
  pattern_branch <- integer()
  for (ci in cells_p) {
    Aci <- fit$A[[ci]]
    for (b in seq_len(R)) {
      nm <- paste0("shift_", labels[b], "_c", ci)
      Ulist[[nm]] <- tcrossprod(Aci[, b])
      pattern_branch[nm] <- b
    }
  }
  ## drift patterns: null covariance of the branch estimates under other phylogenetic parameters
  info <- fit$cell_info[[modal]]
  drift_grid <- if (fit$model == "OUfixedRoot") expand.grid(lambda = drift_lambda, rho_prime = drift_rho_prime) else data.frame(lambda = drift_lambda, rho_prime = 1)
  h <- tree_height(fit$phy_rep)
  for (k in seq_len(nrow(drift_grid))) {
    alpha <- if (fit$model == "OUfixedRoot") rho_prime_inv(drift_grid$rho_prime[k], h, getBoundsSelectionStrength(fit$phy_rep)) else 0
    Ck <- ape::vcv(eb_transform_tree(fit$phy_rep, fit$model, drift_grid$lambda[k], alpha))
    ## excess covariance of the branch estimates relative to the assumed (modal) structure, positive part
    M <- backsolve(info$chol, t(backsolve(info$chol, Ck, transpose = TRUE)), transpose = TRUE) - diag(nrow(Ck))
    U <- crossprod(info$Zp, M %*% info$Zp) / outer(info$nz, info$nz)
    U <- (U + t(U)) / 2
    e <- eigen(U, symmetric = TRUE)
    pos <- e$values > 1e-8 * max(abs(e$values))
    if (!any(pos)) next
    U <- e$vectors[, pos, drop = FALSE] %*% (e$values[pos] * t(e$vectors[, pos, drop = FALSE]))
    Ulist[[paste0("drift_l", drift_grid$lambda[k], if (fit$model == "OUfixedRoot") paste0("_r", drift_grid$rho_prime[k]) else "")]] <- U / max(diag(U))
  }
  if (identity) Ulist <- c(Ulist, mashr::cov_canonical(data, cov_methods = "identity"))
  if (singletons) Ulist <- c(Ulist, mashr::cov_canonical(data, cov_methods = "singletons"))
  if (data_driven) {
    m1 <- mashr::mash_1by1(data)
    strong <- mashr::get_significant_results(m1, lfsr_strong)
    if (length(strong) >= 5 * R) {
      U.pca <- mashr::cov_pca(data, min(n_pca, R), subset = strong)
      U.ed <- mashr::cov_ed(data, U.pca, subset = strong)
      Ulist <- c(Ulist, U.ed)
    } else message("Not enough strong genes for data-driven patterns (", length(strong), "); using canonical patterns only.")
  }
  m <- do.call(mashr::mash, c(list(data = data, Ulist = Ulist, outputlevel = 2, verbose = FALSE), dots))
  lfsr <- ashr::get_lfsr(m); pm <- ashr::get_pm(m); psd <- ashr::get_psd(m)
  ## posterior probability of each single-shift pattern (summed over scales and cells)
  w <- m$posterior_weights
  base <- sub("\\.[0-9]+$", "", colnames(w))
  branch_of <- unname(pattern_branch[base])
  branch_prob <- sapply(seq_len(R), function(b) rowSums(w[, which(branch_of == b), drop = FALSE]))
  if (is.null(dim(branch_prob))) branch_prob <- matrix(branch_prob, nrow = 1)
  dimnames(branch_prob) <- list(fit$genes, labels)
  shift_prob <- rowSums(branch_prob)
  top_branch <- labels[max.col(branch_prob, ties.method = "first")]
  pi <- m$fitted_g$pi
  pattern_weights <- tapply(pi, sub("\\.[0-9]+$", "", names(pi)), sum)
  list(mash = m, lfsr = lfsr, pm = pm, psd = psd, branch_prob = branch_prob,
       shift_prob = shift_prob, top_branch = top_branch,
       pattern_weights = pattern_weights, patterns = names(Ulist), V = Vdat)
}

#' @title Empirical Bayes shift detection on one branch with ash
#'
#' @description
#' Adaptive shrinkage (\code{\link[ashr]{ash}}) of the shift estimates of one
#' branch of the species tree (e.g. the terminal branch of a focal species),
#' from a \code{\link{phyloShiftFit}} object.
#'
#' @param fit a \code{PhyloShiftFit} object.
#' @param branch the label (or index) of the branch to analyse.
#' @param ... further arguments passed to \code{\link[ashr]{ash}}.
#'
#' @return The \code{ash} object; use \code{\link[ashr]{get_lfsr}},
#' \code{\link[ashr]{get_pm}} on it.
#'
#' @export
#'
phyloShiftAsh <- function(fit, branch, ...) {
  if (!requireNamespace("ashr", quietly = TRUE)) stop("Package 'ashr' is needed for function 'phyloShiftAsh'.")
  if (is.character(branch)) branch <- match(branch, fit$branches$label)
  ashr::ash(fit$Bhat[, branch], fit$Shat[, branch], df = fit$df, ...)
}
