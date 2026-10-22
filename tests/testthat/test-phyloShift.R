context("Branch shift detection")

set.seed(7)
S <- 5; reps <- c(2, 3, 4, 5, 6)
tree <- ape::rcoal(S); tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree)); tree$tip.label <- paste0("sp", 1:S)
species <- rep(tree$tip.label, reps)
ids <- paste0(species, "_", ave(seq_along(species), species, FUN = seq_along))
n <- length(ids); G <- 60
Y <- matrix(rnorm(G * n), G, n, dimnames = list(paste0("g", 1:G), ids))
## a clean shift on sp3 for the first 5 genes
Y[1:5, species == "sp3"] <- Y[1:5, species == "sp3"] + 4

test_that("branchIndicators is consistent", {
  bi <- branchIndicators(tree, species)
  expect_equal(dim(bi$Z), c(n, 2 * S - 2))
  expect_equal(sum(bi$branches$testable), 2 * S - 3)
  expect_equal(sum(bi$branches$terminal), S)
  ## terminal indicator = species membership
  expect_equal(unname(bi$Z[, "sp3"]), as.numeric(species == "sp3"))
  ## the root-child branch dropped is the larger clade
  root_edges <- which(bi$branches$parent == S + 1)
  expect_false(bi$branches$testable[root_edges[which.max(bi$branches$n_species[root_edges])]])
  ## internal node indicator = union of its children
  int <- which(!bi$branches$terminal & bi$branches$testable)[1]
  kids <- which(bi$branches$parent == bi$branches$child[int])
  expect_equal(bi$Z[, int], rowSums(bi$Z[, kids, drop = FALSE]))
})

test_that("phyloShiftFit recovers a clean shift and matches phylolmFit", {
  sf <- phyloShiftFit(Y, phy = tree, col_species = species, model = "BM", n_lambda = 9)
  expect_s3_class(sf, "PhyloShiftFit")
  expect_equal(dim(sf$Bhat), c(G, 2 * S - 3))
  expect_true(all(sf$p.value[1:5, "sp3"] < 1e-4))
  expect_true(mean(sf$p.value[6:G, ] < 0.05) < 0.2)
  expect_equal(sf$t, sf$Bhat / sf$Shat)
  ## plug-in variant equals a phylolmFit with the indicator in the design at the same grid cell
  sfp <- phyloShiftFit(Y, phy = tree, col_species = species, model = "BM", n_lambda = 9, nuisance = "eb_plugin", squeeze = FALSE)
  g <- 3; b <- "sp3"
  grid <- sfp$eb_null$grid
  cell <- snap_to_grid(sfp$lambda[g, b], 1, grid)
  X <- cbind(1, sfp$Z[, b]); dimnames(X) <- list(ids, c("(Intercept)", "branch"))
  eb1 <- sfp$eb_null; eb1$cell <- cell; eb1$genes <- rownames(Y)[g]
  f1 <- phylolmFit(Y[g, , drop = FALSE], design = X, phy = sfp$phy_rep, model = "BM", eb_model = eb1)
  expect_equal(unname(sfp$Bhat[g, b]), unname(f1$coefficients[, "branch"]), tolerance = 1e-8)
  expect_equal(unname(sfp$Shat[g, b]), unname(f1$stdev.unscaled[, "branch"] * f1$sigma), tolerance = 1e-8)
  ## leakage: a unit shift on a terminal branch is diluted in its parent's test
  A <- sf$A[[1]]
  expect_equal(unname(diag(A)), rep(1, ncol(A)))
  ## ash on one branch
  a <- phyloShiftAsh(sf, "sp3")
  expect_true(all(ashr::get_lfsr(a)[1:5] < 0.01))
})

test_that("phyloMash recovers a clean shift with calibrated probabilities", {
  expect_warning(pm <- phyloMash(Y, tree, col_species = species, model = "BM", n_grid = 8), "fitted on")
  expect_s3_class(pm, "PhyloMash")
  expect_equal(dim(pm$branch_prob), c(G, 2 * S - 3))
  expect_true(all(pm$branch_prob[1:5, "sp3"] > 0.9))
  expect_true(all(pm$top_branch[1:5] == "sp3"))
  expect_true(all(abs(rowSums(pm$pm_species)) < 1e-8))
  calls <- callShifts(pm, 0.05)
  expect_true(all(calls[1:5, "sp3"]))
  expect_true(sum(calls[6:G, ]) <= 3)
  expect_equal(sum(pm$pattern_weights), 1, tolerance = 1e-6)
  ## OU drift patterns and subset fitting run
  expect_warning(pm2 <- phyloMash(Y, tree, col_species = species, model = "OUfixedRoot", n_grid = 6, fit_subset = 40), "fitted on")
  expect_true(any(grepl("^OU", names(pm2$prior_fit$drift_shape))))
  pm3 <- phyloMash(Y, tree, col_species = species, model = "OUfixedRoot", prior = "mash", n_grid = 6, fit_subset = 40)
  expect_true(any(grepl("drift_OU", names(pm3$pattern_weights))))
})

test_that("phyloMash shift probabilities count shifts above min_fold", {
  expect_warning(pm1 <- phyloMash(Y, tree, col_species = species, model = "BM", n_grid = 8, min_fold = 1), "fitted on")
  expect_warning(pm2 <- phyloMash(Y, tree, col_species = species, model = "BM", n_grid = 8, min_fold = 1.25), "fitted on")
  ## the threshold changes the probabilities, not the fitted prior
  expect_equal(pm1$fit$pi, pm2$fit$pi)
  expect_equal(factorized_posterior(pm1$data, pm1$fit, min_shift = log2(1.25))$branch_prob, unname(pm2$branch_prob))
  expect_true(all(pm2$branch_prob <= pm1$branch_prob + 1e-12))
  expect_true(all(pm2$branch_prob[1:5, "sp3"] > 0.9))                              # shift of 4 (log2) on sp3
  expect_lt(max(factorized_posterior(pm1$data, pm1$fit, min_shift = 10)$branch_prob), 1e-3)
  expect_error(phyloMash(Y, tree, col_species = species, min_fold = 0.8), "min_fold")
  ## re-scoring without refitting gives the fit at that threshold
  pm3 <- rescoreShifts(pm1, 1.25)
  expect_equal(pm3$branch_prob, pm2$branch_prob); expect_equal(pm3$shift_prob, pm2$shift_prob)
  expect_equal(pm3$top_branch, pm2$top_branch); expect_equal(pm3$min_fold, 1.25)
  expect_equal(rescoreShifts(pm3, 1)$branch_prob, pm1$branch_prob)
  expect_error(rescoreShifts(pm1, 0.5), "min_fold")
  ## the normal posterior of the shift size given a component matches direct conditioning
  f <- pm1$fit; d <- pm1$data
  ci <- which(f$comp$b > 0 & f$pi > 1e-6)[1]; a <- f$dirs[[f$comp$b[ci]]]; v <- f$v_grid[f$comp$v[ci]]
  qf <- shift_quadforms(d$Bhat, d$Shat, d$V, f$U[[(ci - 1) %% f$K3 + 1]], cbind(a))
  for (g in c(1, 7, 30)) {
    Sig <- diag(d$Shat[g, ]) %*% d$V %*% diag(d$Shat[g, ]) + f$U[[ci]]
    expect_equal(v * qf$r[g, 1] / (1 + v * qf$p[g, 1]), drop(v * t(a) %*% solve(Sig, d$Bhat[g, ])), tolerance = 1e-8)
    expect_equal(v / (1 + v * qf$p[g, 1]), drop(v - v^2 * t(a) %*% solve(Sig, a)), tolerance = 1e-8)
  }
})
