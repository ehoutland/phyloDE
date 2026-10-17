context("Weighted phylogenetic fit and cross analysis")

set.seed(3)
S <- 6; reps <- c(3, 3, 4, 8, 4, 4)
tree <- ape::rcoal(S); tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree)); tree$tip.label <- paste0("sp", 1:S)
species <- rep(tree$tip.label, reps)
case <- as.numeric(species == "sp4" & ave(seq_along(species), species, FUN = seq_along) > 4)   # last 4 of sp4's 8 are cases
ids <- paste0(species, ifelse(case == 1, "_case", ""), "_", ave(seq_along(species), species, case, FUN = seq_along))
n <- length(ids); G <- 120
mu <- rnorm(G, log(40), 1); lib <- exp(rnorm(n, 0, 1))
counts <- matrix(rnbinom(G * n, mu = outer(exp(mu), lib), size = 4), G, n, dimnames = list(paste0("g", 1:G), ids))
counts[1:6, species == "sp2"] <- rnbinom(6 * sum(species == "sp2"), mu = outer(exp(mu[1:6]), lib[species == "sp2"]) * 8, size = 4)  # shift on sp2
counts[7:12, case == 1] <- rnbinom(6 * sum(case == 1), mu = outer(exp(mu[7:12]), lib[case == 1]) * 8, size = 4)                 # case effect

test_that("voomPhyloWeights tracks depth and phyloWeightedFit recovers a shift", {
  v <- voomPhyloWeights(counts, species)
  expect_equal(dim(v$weights), dim(counts))
  expect_true(all(v$weights > 0))
  expect_gt(median(apply(log(v$weights), 1, function(w) cor(w, log(lib)))), 0.5)
  wf <- phyloWeightedFit(counts, tree, species, model = "BM", n_lambda = 9)
  expect_s3_class(wf, "PhyloShiftFit")
  expect_true(all(wf$p.value[1:6, "sp2"] < 1e-3))
  expect_true(mean(wf$p.value[13:G, ] < 0.05) < 0.15)
  ## unit weights give the unweighted (sample-level) model
  uf <- phyloWeightedFit(counts, tree, species, model = "BM", n_lambda = 9, weights = matrix(1, G, n, dimnames = dimnames(counts)))
  expect_equal(dim(uf$Bhat), dim(wf$Bhat))
  ## exact vs profile-shared factors agree
  ef <- phyloWeightedFit(counts, tree, species, model = "BM", n_lambda = 9, n_profiles = Inf)
  expect_gt(cor(as.vector(wf$t), as.vector(ef$t)), 0.99)
  ## design coefficient
  d <- cbind(1, case); dimnames(d) <- list(ids, c("(Intercept)", "case"))
  cf <- phyloWeightedFit(counts, tree, species, design = d, model = "BM", n_lambda = 9, shifts = FALSE)
  expect_equal(colnames(cf$Bhat), "case")
  expect_true(all(cf$p.value[7:12, "case"] < 1e-3))
})

test_that("speciesStatsWeighted feeds phyloMash", {
  st <- speciesStatsWeighted(counts, species)
  expect_s3_class(st, "SpeciesStats"); expect_equal(st$method, "weighted")
  expect_true(all(is.finite(st$se)))
  pm <- phyloMash(st, tree, model = "BM", n_grid = 8)
  expect_true(all(pm$branch_prob[1:6, "sp2"] > 0.9))
  pc <- phyloMash(st, tree, model = "BM", n_grid = 8, patterns = "canonical")
  expect_equal(pc$patterns, "canonical")
  expect_equal(dim(pc$branch_prob), dim(pm$branch_prob))
})

test_that("phyloMashCross classifies evo-only and case-only genes", {
  cr <- phyloMashCross(counts, tree, species, case, focal = "sp4", model = "BM", n_lambda = 9)
  expect_s3_class(cr, "PhyloMashCross")
  expect_equal(colnames(cr$Bhat), c("evo", "case"))
  expect_true(all(abs(rowSums(cr$class_prob) - 1) < 1e-8))
  expect_true(all(cr$top_class[7:12] == "case_only"))
  expect_true(all(cr$class_prob[7:12, "case_only"] > 0.8))
  expect_true(mean(cr$class_prob[13:G, "null"]) > mean(cr$class_prob[1:12, "null"]))
  cl <- classifyCross(cr, 0.05)
  expect_equal(colnames(cl), c("evo_only", "case_only", "shared_same", "shared_opposite", "independent"))
  expect_true(all(cl[7:12, "case_only"]))
  expect_true(mean(rowSums(cl[13:G, ]) > 0) < 0.15)   # few null genes called in any class
  expect_true(abs(cr$V[1, 2]) < 1)
  expect_error(phyloMashCross(counts, tree, species, as.numeric(species == "sp1"), focal = "sp4"), "focal")
})
