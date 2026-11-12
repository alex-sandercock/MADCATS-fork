#!/usr/bin/env Rscript
# Writes tests/fixtures/polyrad_parity.json: reference values from polyRAD that
# tests/test_call.py compares the madcats_call genotype engine against.
#
# Usage (from the repository root, with polyRAD installed):
#   Rscript scripts/make_polyrad_fixture.R
#
# The fixture holds two sets of values for one small simulated data set:
#   single_pass  one round of likelihood, HWE prior and posterior, starting
#                from the depth-ratio allele frequencies.  madcats_call must
#                reproduce these almost exactly.
#   iterated     IterateHWE run to a very tight tolerance, plus polyRAD's most
#                probable genotypes after its correction to the ploidy.
# The committed JSON was made with the polyRAD version recorded inside it.

suppressMessages(library(polyRAD))

overdispersion <- 9
contamination_rate <- 0.001
iteration_tolerance <- 1e-10

# Three loci with 2, 3 and 4 alleles; 16 samples of ploidy 2, 3, 4 and 6.
set.seed(20240611)
alleles_per_locus <- c(2L, 3L, 4L)
alleles2loc <- rep(seq_along(alleles_per_locus), alleles_per_locus)
sample_ploidy <- rep(c(2L, 3L, 4L, 6L), each = 4)
sample_names <- sprintf("sample%02d", seq_along(sample_ploidy))
allele_names <- unlist(lapply(seq_along(alleles_per_locus), function(locus) {
  sprintf("locus%d_allele%d", locus, seq_len(alleles_per_locus[locus]))
}))

depth <- matrix(0L, nrow = length(sample_names), ncol = length(allele_names),
                dimnames = list(sample_names, allele_names))
for (locus in seq_along(alleles_per_locus)) {
  columns <- which(alleles2loc == locus)
  frequencies <- rev(seq_along(columns)) / sum(seq_along(columns))
  for (sample in seq_along(sample_names)) {
    copies <- rmultinom(1, sample_ploidy[sample], frequencies)[, 1]
    total_reads <- sample(c(4:60, 150, 400), 1)
    depth[sample, columns] <- rmultinom(1, total_reads, copies + 0.02)[, 1]
  }
}
# One sample with no reads at the first locus, and one with a single read.
depth["sample03", alleles2loc == 1] <- 0L
depth["sample09", alleles2loc == 1] <- c(1L, 0L)

loc_table <- data.frame(row.names = sprintf("locus%d", seq_along(alleles_per_locus)))
nucleotides <- unlist(lapply(alleles_per_locus, function(n) c("A", "C", "G", "T")[seq_len(n)]))
rad <- RADdata(depth, alleles2loc, loc_table, list(2L), contamination_rate,
               nucleotides, setNames(sample_ploidy, sample_names))

# Flatten polyRAD's [dosage, sample, allele] arrays (one per sample ploidy)
# into one record per sample and allele.
records <- function(object, with_likelihood) {
  out <- list()
  for (h in seq_len(ncol(object$posteriorProb))) {
    posterior <- object$posteriorProb[[1, h]]
    for (sample in dimnames(posterior)[[2]]) {
      for (allele in dimnames(posterior)[[3]]) {
        record <- list(sample = sample, allele = allele,
                       posterior = unname(posterior[, sample, allele]))
        if (with_likelihood) {
          record$likelihood <- unname(object$genotypeLikelihood[[1, h]][, sample, allele])
          record$prior <- unname(object$priorProb[[1, h]][, allele])
        }
        out[[length(out) + 1]] <- record
      }
    }
  }
  out
}

single <- AddAlleleFreqHWE(rad)
start_frequency <- single$alleleFreq
single <- AddGenotypePriorProb_HWE(single)
single <- AddGenotypeLikelihood(single, overdispersion = overdispersion)
single <- AddPloidyChiSq(single)
single <- AddGenotypePosteriorProb(single)

iterated <- suppressMessages(IterateHWE(rad, tol = iteration_tolerance,
                                        overdispersion = overdispersion))
genotypes <- GetProbableGenotypes(iterated, omit1allelePerLocus = FALSE,
                                  multiallelic = "correct")$genotypes

# Minimal JSON writer, so the script needs no package beyond polyRAD.
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
  polyrad_version = as.character(packageVersion("polyRAD")),
  overdispersion = overdispersion,
  contamination_rate = contamination_rate,
  iteration_tolerance = iteration_tolerance,
  samples = as.list(sample_names),
  ploidy = as.list(sample_ploidy),
  alleles = as.list(allele_names),
  allele_locus = as.list(alleles2loc),
  depth = lapply(seq_len(nrow(depth)), function(row) unname(depth[row, ])),
  single_pass = list(
    allele_frequency = unname(start_frequency),
    records = records(single, TRUE)
  ),
  iterated = list(
    allele_frequency = unname(iterated$alleleFreq),
    records = records(iterated, FALSE),
    genotypes = lapply(seq_len(nrow(genotypes)), function(row) unname(genotypes[row, ]))
  )
)
writeLines(to_json(fixture), "tests/fixtures/polyrad_parity.json")
cat("Wrote tests/fixtures/polyrad_parity.json\n")
