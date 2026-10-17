context("Species-level statistics and GLS")

set.seed(11)
S <- 5; reps <- c(3, 3, 4, 5, 6)
tree <- ape::rcoal(S); tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree)); tree$tip.label <- paste0("sp", 1:S)
species <- rep(tree$tip.label, reps)
ids <- paste0(species, "_", ave(seq_along(species), species, FUN = seq_along))
n <- length(ids); G <- 80
mu <- rnorm(G, log(50), 1)
counts <- matrix(rnbinom(G * n, mu = exp(mu), size = 5), G, n, dimnames = list(paste0("g", 1:G), ids))
counts[1:5, species == "sp3"] <- rnbinom(5 * sum(species == "sp3"), mu = exp(mu[1:5]) * 8, size = 5)   # shift on sp3
counts[6, species == "sp1"] <- 0                                                                        # dropout
Y <- lengthNormalizeRNASeq(counts, lengthNormalization = "none", dataTransformation = "log2")

test_that("speciesStats limma and edgeR agree and flag dropouts", {
  sl <- speciesStats(Y, species, method = "limma")
  se <- speciesStats(counts, species, method = "edgeR")
  expect_s3_class(sl, "SpeciesStats"); expect_s3_class(se, "SpeciesStats")
  expect_equal(dim(sl$mean), c(G, S)); expect_equal(dim(se$se), c(G, S))
  expect_true(all(is.finite(se$se))); expect_true(all(is.finite(se$mean)))
  expect_gt(cor(as.vector(sl$mean[-6, ]), as.vector(se$mean[-6, ])), 0.98)
  expect_true(se$dropout[6, "sp1"]); expect_false(any(se$dropout[-6, ]))
  ## more replicates -> smaller SE, both methods
  expect_true(median(se$se[, "sp5"]) < median(se$se[, "sp1"]))
  expect_true(median(sl$se[, "sp5"]) < median(sl$se[, "sp1"]))
  ## species means equal the group means of Y
  expect_equal(unname(sl$mean[3, "sp2"]), mean(Y[3, species == "sp2"]))
  ## covariates: projecting out a covariate orthogonal to species leaves the means unchanged
  cov1 <- rnorm(n); cov1 <- cov1 - ave(cov1, species)
  slc <- speciesStats(Y, species, method = "limma", covariates = cbind(cov1))
  expect_equal(slc$mean, sl$mean, tolerance = 1e-10)
})

test_that("phyloShiftGLS recovers the shift and is calibrated", {
  se <- speciesStats(counts, species, method = "edgeR")
  design <- cbind(1, as.numeric(tree$tip.label %in% c("sp3", "sp4"))); dimnames(design) <- list(tree$tip.label, c("(Intercept)", "cond"))
  g <- phyloShiftGLS(se, tree, design = design, model = "BM", n_lambda = 9)
  expect_s3_class(g, "PhyloShiftFit")
  expect_equal(colnames(g$Bhat)[1], "cond")
  expect_true(all(g$p.value[1:5, "sp3"] < 1e-3))
  expect_true(all(abs(g$Bhat[1:5, "sp3"] - 3) < 1.2))   # log2(8) = 3
  expect_true(mean(g$p.value[7:G, ] < 0.05) < 0.15)
  ## OU runs
  g2 <- phyloShiftGLS(se, tree, model = "OUfixedRoot", n_lambda = 7, n_rho = 4, shifts = FALSE, design = design)
  expect_equal(ncol(g2$Bhat), 1)
})

test_that("phyloMash accepts SpeciesStats input", {
  se <- speciesStats(counts, species, method = "edgeR")
  pm <- phyloMash(se, tree, model = "BM", n_grid = 8)
  expect_s3_class(pm, "PhyloMash")
  expect_equal(pm$stats_method, "edgeR")
  expect_true(all(pm$branch_prob[1:5, "sp3"] > 0.9))
  expect_equal(dim(pm$dropout), c(G, S))
  pm2 <- phyloMash(Y, tree, col_species = species, model = "BM", n_grid = 8)
  expect_equal(pm2$stats_method, "limma")
  expect_true(all(pm2$branch_prob[1:5, "sp3"] > 0.9))
})
