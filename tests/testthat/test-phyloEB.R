context("Empirical Bayes grid estimation of phylogenetic parameters")

data(crayfish)
norm_data <- lengthNormalizeRNASeq(crayfish$counts, crayfish$lengths)
design <- model.matrix(~ sights, model.frame(crayfish$sights))
y <- check_expression_matrix(norm_data[1:20, ])
phy <- check_tree(crayfish$tree, y, NULL)
y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
design <- check_design_matrix(design, y, phy)

test_that("grid REML log-likelihood matches phylolm at fixed lambda (BM + ME)", {
  grid <- phyloEBGrid(phy, "BM", TRUE, n_lambda = 7)
  LL <- phyloEBLogLik(y_data, design, phy, "BM", grid, REML = TRUE)
  LLml <- phyloEBLogLik(y_data, design, phy, "BM", grid, REML = FALSE)
  for (k in 2:6) {
    for (g in 1:3) {
      dd <- data.frame(expr = y_data[g, ], design)
      fit <- phylolm::phylolm(expr ~ -1 + ., data = dd, phy = phy, model = "lambda", REML = TRUE,
                              lower.bound = list(lambda = grid$lambda[k]), upper.bound = list(lambda = grid$lambda[k]),
                              starting.value = list(lambda = grid$lambda[k]))
      expect_equal(unname(LL[g, k]), fit$logLik, tolerance = 1e-8)
      fitml <- phylolm::phylolm(expr ~ -1 + ., data = dd, phy = phy, model = "lambda", REML = FALSE,
                                lower.bound = list(lambda = grid$lambda[k]), upper.bound = list(lambda = grid$lambda[k]),
                                starting.value = list(lambda = grid$lambda[k]))
      expect_equal(unname(LLml[g, k]), fitml$logLik, tolerance = 1e-8)
    }
  }
})

test_that("grid REML log-likelihood matches phylolm at its optimum (OU + ME)", {
  alpha <- 2 / tree_height(phy)
  tree_ou <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha))$tree
  tilde_t <- tree_height(tree_ou) / (2 * alpha)
  for (g in 1:3) {
    dd <- data.frame(expr = y_data[g, ], design)
    fit <- suppressWarnings(phylolm::phylolm(expr ~ -1 + ., data = dd, phy = phy, model = "OUfixedRoot", measurement_error = TRUE, REML = TRUE,
                                             lower.bound = list(alpha = alpha), upper.bound = list(alpha = alpha), starting.value = list(alpha = alpha)))
    lambda_hat <- fit$sigma2 * tilde_t / (fit$sigma2_error + fit$sigma2 * tilde_t)
    grid <- data.frame(lambda = lambda_hat, rho_prime = rho_prime(alpha, tree_height(phy)), alpha = alpha)
    LL <- phyloEBLogLik(y_data[g, , drop = FALSE], design, phy, "OUfixedRoot", grid, REML = TRUE)
    expect_equal(as.vector(LL), fit$logLik, tolerance = 1e-6)
  }
})

test_that("phyloEB and phylolmFit with eb_model run end to end", {
  eb <- phyloEB(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "BM", n_lambda = 11)
  expect_s4_class(eb, "PhyloEBModel")
  expect_equal(sum(eb$prior), 1)
  expect_true(all(abs(rowSums(eb$posterior) - 1) < 1e-10))
  expect_equal(length(eb$lambda), 20)
  expect_true(all(eb$lambda >= 0 & eb$lambda <= 1))
  pp <- getParameters(eb, consensus = FALSE)
  expect_equal(colnames(pp), "lambda")
  expect_equal(names(getParameters(eb)), "lambda")

  fit <- phylolmFit(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "BM", eb_model = eb)
  expect_s4_class(fit, "PhyloMArrayLM")
  expect_false(fit$use_consensus)
  expect_equal(fit$lambda_error, eb$grid$lambda[eb$cell])
  fitb <- limma::eBayes(fit)
  tt <- limma::topTable(fitb, coef = 2, number = Inf)
  expect_equal(nrow(tt), 20)
  expect_equal(getParameters(fit, consensus = FALSE)$lambda, eb$grid$lambda[eb$cell])

  ## the EB fit for a gene equals a gene-wise fit at the snapped lambda
  g <- 3
  tr <- eb_transform_tree(crayfish$tree, "BM", eb$grid$lambda[eb$cell[g]])
  ## replicate tree is needed
  tr <- eb_transform_tree(phy, "BM", eb$grid$lambda[eb$cell[g]])
  R <- chol(ape::vcv(tr))
  yt <- backsolve(R, y_data[g, ], transpose = TRUE)
  Xt <- backsolve(R, design, transpose = TRUE)
  ols <- lm.fit(Xt, yt)
  expect_equal(unname(fit$coefficients[g, ]), unname(ols$coefficients), tolerance = 1e-8)

  ## OU
  ebou <- phyloEB(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "OUfixedRoot", n_lambda = 7, n_rho = 5)
  expect_equal(colnames(getParameters(ebou, consensus = FALSE)), c("lambda", "rho"))
  fitou <- phylolmFit(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "OUfixedRoot", eb_model = ebou)
  expect_equal(fitou$optpar, ebou$grid$alpha[ebou$cell])
  expect_silent(limma::eBayes(fitou))

  ## errors
  expect_error(phylolmFit(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "OUfixedRoot", eb_model = eb), "model")
  expect_error(phylolmFit(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "BM", eb_model = eb, use_consensus = TRUE), "use_consensus")
  expect_error(phylolmFit(norm_data[2:21, ], design = design, phy = crayfish$tree, model = "BM", eb_model = eb), "genes")
  expect_error(phyloEB(norm_data[1:20, ], design = design, phy = crayfish$tree, model = "BM", measurement_error = FALSE), "no phylogenetic parameter")
})

test_that("prior smoothing components are proper", {
  grid <- phyloEBGrid(phy, "OUfixedRoot", TRUE, n_lambda = 5, n_rho = 4)
  B <- eb_prior_components(grid, c(0.5, 2))
  expect_equal(dim(B), c(20, 40))
  expect_true(all(abs(colSums(B) - 1) < 1e-12))
  expect_equal(eb_prior_components(grid, NULL), diag(20))
})
