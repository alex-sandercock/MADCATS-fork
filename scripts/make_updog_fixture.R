#!/usr/bin/env Rscript
# Writes tests/fixtures/updog_norm.json: reference values from updog that
# tests/test_call.py compares the normal-shaped prior of madcats_call against.
#
# Usage (from the repository root, with updog installed):
#   Rscript scripts/make_updog_fixture.R
#
# The fixture holds two sets of values:
#   prior_fits  updog's fit of the norm genotype distribution to given dosage
#               weights (its internal flex_update_pivec), for weights of many
#               shapes.  madcats_call.fit_normal_prior must find the same
#               distribution.
#   markers     complete updog fits (flexdog_full, model = "norm") of simulated
#               read counts, with updog's allelic bias, sequencing error and
#               overdispersion held fixed at no bias, no error and the given
#               overdispersion.  With those fixed, updog's likelihood is the
#               beta-binomial that madcats_call uses with a contamination rate
#               of 0, so the posteriors must agree.
# The committed JSON was made with the updog version recorded inside it.

suppressMessages(library(updog))

# madcats_call and polyRAD give the overdispersion as the sum of the two
# beta-binomial shape parameters; updog as rho = 1 / (that sum + 1).
overdispersion <- 30
updog_overdispersion <- 1 / (overdispersion + 1)

weight_sets <- list(
  list(ploidy = 4L, weights = c(10, 40, 30, 15, 5)),            # one peak inside the range
  list(ploidy = 4L, weights = c(90, 9, 0.9, 0.09, 0.01)),       # rare allele: mean below 0
  list(ploidy = 4L, weights = c(0.01, 0.5, 3, 30, 66)),         # common allele
  list(ploidy = 4L, weights = c(50, 0.5, 0.2, 0.5, 50)),        # two peaks: inbred lines
  list(ploidy = 4L, weights = c(20, 20, 20, 20, 20)),           # flat
  list(ploidy = 4L, weights = c(0, 1, 98, 1, 0)),               # one F1 class
  list(ploidy = 4L, weights = c(60, 40, 0.001, 0.001, 0.001)),  # two neighbouring dosages
  list(ploidy = 2L, weights = c(30, 50, 20)),
  list(ploidy = 2L, weights = c(49, 2, 49)),
  list(ploidy = 2L, weights = c(77, 51, 45)),
  list(ploidy = 6L, weights = c(1, 2, 5, 20, 40, 25, 7)),
  list(ploidy = 6L, weights = c(61, 13, 7, 21, 8, 0.4, 2.8)),
  list(ploidy = 3L, weights = c(5, 30, 40, 10)),
  list(ploidy = 8L, weights = c(2, 9, 20, 30, 22, 11, 4, 1.5, 0.5))
)
prior_fits <- lapply(weight_sets, function(set) {
  fit <- updog:::flex_update_pivec(weight_vec = set$weights, model = "norm",
                                   control = list(ploidy = set$ploidy, mu = set$ploidy / 2, sigma = 1))
  list(ploidy = set$ploidy, weights = set$weights, prior = as.vector(fit$pivec),
       mean = fit$par$mu, standard_deviation = fit$par$sigma)
})

# Simulated read counts of one allele: genotype classes of three kinds of population.
set.seed(20240930)
simulate <- function(ploidy, dosages) {
  depth <- sample(15:90, length(dosages), replace = TRUE)
  fraction <- dosages / ploidy
  shape1 <- pmax(fraction * overdispersion, 1e-9); shape2 <- pmax((1 - fraction) * overdispersion, 1e-9)
  read_fraction <- ifelse(fraction %in% c(0, 1), fraction, rbeta(length(dosages), shape1, shape2))
  list(ploidy = ploidy, allele_depth = rbinom(length(dosages), depth, read_fraction), marker_depth = depth)
}
marker_sets <- list(
  random_mating_tetraploid = simulate(4L, rbinom(90, 4, 0.3)),
  one_family_hexaploid     = simulate(6L, sample(c(2, 3, 4), 80, replace = TRUE, prob = c(0.25, 0.5, 0.25))),
  inbred_diploid           = simulate(2L, sample(c(0, 2), 70, replace = TRUE, prob = c(0.6, 0.4))),
  rare_allele_tetraploid   = simulate(4L, rbinom(90, 4, 0.03))
)
markers <- lapply(names(marker_sets), function(name) {
  data <- marker_sets[[name]]
  fit <- flexdog_full(refvec = data$allele_depth, sizevec = data$marker_depth, ploidy = data$ploidy,
                      model = "norm", verbose = FALSE, seq = 0, bias = 1, od = updog_overdispersion,
                      update_bias = FALSE, update_seq = FALSE, update_od = FALSE,
                      itermax = 20000, tol = 1e-13)
  list(name = name, ploidy = data$ploidy, allele_depth = data$allele_depth, marker_depth = data$marker_depth,
       prior = as.vector(fit$gene_dist), mean = fit$par$mu, standard_deviation = fit$par$sigma,
       genotypes = as.vector(fit$geno),
       posteriors = lapply(seq_len(nrow(fit$postmat)), function(row) as.vector(fit$postmat[row, ])))
})

# Minimal JSON writer, so the script needs no package beyond updog.
to_json <- function(value) {
  if (is.list(value)) {
    if (is.null(names(value))) {
      return(paste0("[", paste(vapply(value, to_json, ""), collapse = ","), "]"))
    }
    fields <- vapply(names(value), function(name) {
      paste0("\"", name, "\":", to_json(value[[name]]))
    }, "")
    return(paste0("{", paste(fields, collapse = ","), "}"))
  }
  if (is.character(value)) {
    text <- paste0("\"", value, "\"")
  } else {
    text <- ifelse(is.na(value), "null", sprintf("%.17g", value))
  }
  if (length(value) == 1) text else paste0("[", paste(text, collapse = ","), "]")
}

fixture <- list(
  updog_version = as.character(packageVersion("updog")),
  overdispersion = overdispersion,
  prior_fits = prior_fits,
  markers = markers
)
writeLines(to_json(fixture), "tests/fixtures/updog_norm.json")
cat("Wrote tests/fixtures/updog_norm.json\n")
