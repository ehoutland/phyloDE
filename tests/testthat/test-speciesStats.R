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

test_that("drift_prior = 0 is the plain mash fit and the default lowers the shift weight", {
  se <- speciesStats(counts, species, method = "edgeR")
  pm0 <- phyloMash(se, tree, model = "BM", prior = "mash", n_grid = 8, drift_prior = 0)
  g0 <- pm0$mash$fitted_g
  data <- mashr::mash_set_data(pm0$species_means %*% t(pm0$L), sqrt(pm0$se_species^2 %*% t(pm0$L^2)), df = pm0$df,
                               V = stats::cov2cor(pm0$L %*% diag(colMeans(pm0$se_species^2) / mean(pm0$se_species^2)) %*% t(pm0$L)))
  m_ref <- mashr::mash(data, Ulist = g0$Ulist, grid = g0$grid, outputlevel = 1, verbose = FALSE)
  expect_equal(unname(g0$pi), unname(m_ref$fitted_g$pi), tolerance = 1e-6)
  shift_w <- function(pm) sum(pm$pattern_weights[grepl("^shift_", names(pm$pattern_weights)) & !grepl("_r1$", names(pm$pattern_weights))])
  pm1 <- phyloMash(se, tree, model = "BM", prior = "mash", n_grid = 8, drift_prior = 0.5)
  expect_lt(shift_w(pm1), shift_w(pm0))
  expect_equal(pm1$drift_prior, 0.5)
  expect_true(all(pm1$branch_prob[1:5, "sp3"] > 0.9))   # a strong shift survives a heavy penalty
})

test_that("sample quality weights downweight outlier libraries in both routes", {
  set.seed(4)
  C <- counts; bad <- c(2, 9)   # one library of sp1 and one of sp3 with extra noise
  C[, bad] <- round(C[, bad] * exp(matrix(rnorm(G * 2, 0, 1), ncol = 2)))
  w <- sampleQualityWeights(C, species)
  expect_length(w, n); expect_named(w, ids)
  expect_true(all(w[bad] < 0.5) && median(w[-bad]) > 0.8)
  for (m in c("limma", "edgeR")) {
    norm <- function(x) if (m == "limma") lengthNormalizeRNASeq(x, lengthNormalization = "none", dataTransformation = "log2") else x
    ref <- speciesStats(norm(counts), species, method = m)   # clean data, same route
    s0 <- speciesStats(norm(C), species, method = m); s1 <- speciesStats(norm(C), species, method = m, sample_weights = w)
    expect_equal(dim(s1$se), dim(s0$se))
    ## weighting brings the means of the affected species closer to the clean-data means
    for (s in c("sp1", "sp3")) expect_lt(mean((s1$mean[, s] - ref$mean[, s])^2), mean((s0$mean[, s] - ref$mean[, s])^2))
  }
  expect_error(speciesStats(C, species, method = "edgeR", sample_weights = w[-1]), "sample_weights")
  wf <- phyloWeightedFit(C, tree, species, model = "BM", n_lambda = 7, sample_weights = w)
  expect_s3_class(wf, "PhyloShiftFit")
})
