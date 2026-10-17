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
  pm <- phyloMash(Y, tree, col_species = species, model = "BM", n_grid = 8)
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
  pm2 <- phyloMash(Y, tree, col_species = species, model = "OUfixedRoot", n_grid = 6, fit_subset = 40)
  expect_true(any(grepl("drift_OU", names(pm2$pattern_weights))))
})
