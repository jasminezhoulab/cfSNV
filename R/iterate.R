

collect_indel_related_regions <- function(sample.id, reference, python) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  rebase_dir <- paste0(tmp.dir, "rebase/")

  estimate_dir <- paste0(tmp.dir, "estimate/")
  if (system.file("extdata/tmp/estimate", package = "cfSNV") == "") {
    system2(command = "mkdir", args = estimate_dir)
  }

  variant_dir <- paste0(tmp.dir, "variant/")
  if (system.file("extdata/tmp/variant", package = "cfSNV") == "") {
    system2(command = "mkdir", args = variant_dir)
  }

  filter_dir <- paste0(tmp.dir, "filter/")
  if (system.file("extdata/tmp/filter", package = "cfSNV") == "") {
    system2(command = "mkdir", args = filter_dir)
  }

  ## Collect regions affected by indels ---------------------------------------------------
  indel <- paste0(rebase_dir, sample.id, ".indel")
  filter_indel.bed <- paste0(rebase_dir, sample.id, ".filter_indel.bed")

  py14 <- paste0(python.script.dir, "/py14.collect_indel_related_regions.py")
  py14.command <- paste(py14, indel, filter_indel.bed, paste0(reference, ".fai"))
  system2(command = python, args = py14.command)


}


iterate_filter <- function(SNP.database, depth, threshold, sample.id, python, MIN_PASS_SUPPORT_COUNT) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  rebase_dir <- paste0(tmp.dir, "rebase/")
  estimate_dir <- paste0(tmp.dir, "estimate/")
  variant_dir <- paste0(tmp.dir, "variant/")
  filter_dir <- paste0(tmp.dir, "filter/")

  rebase <- paste0(rebase_dir, sample.id, ".rebase")

  ## MUTATION CLUSTER FREQUENCY ESTIMATION -----------------------------------------------
  C_merged_normal.HOTSPOT.VAF <- paste0(
    estimate_dir, sample.id, ".", threshold, ".C_merged_normal.HOTSPOT.VAF")
  C_merged_normal.estimate <- paste0(
    estimate_dir, sample.id, ".", threshold, ".C_merged_normal.estimate")

  # py2 <- paste0(python.script.dir, "/py2.estimate.TFestimate.py")
  # py2.command <- paste(py2, rebase, threshold, sample.id, depth,
  #                      C_merged_normal.HOTSPOT.VAF, C_merged_normal.estimate)
  # system2(command = python, args = py2.command)

  estimate_TFestimate_main(rebase, threshold, sample.id, depth, C_merged_normal.HOTSPOT.VAF, C_merged_normal.estimate)
  estimate <- as.numeric(readLines(C_merged_normal.estimate, n = 1))

  ## CANDIDATE DETECTION AT CURRENT FREQUENCY --------------------------------------------
  variant <- paste0(variant_dir, sample.id, ".", threshold, ".variant")

  # py3 <- paste0(python.script.dir, "/py3.genotype.genotype.py")
  # py3.command <- paste(py3, rebase, variant, estimate, threshold, depth)
  # system2(command = python, args = py3.command)

  genotype_genotype_main(rebase, variant, estimate, threshold, depth, MIN_PASS_SUPPORT_COUNT)

  system2(command = "rm", args = C_merged_normal.HOTSPOT.VAF)
  system2(command = "rm", args = C_merged_normal.estimate)

  ## SITE-LEVEL FILTRATION ---------------------------------------------------------------
  intermediate <- paste0(filter_dir, sample.id, ".", threshold, ".intermediate")
  pass <- paste0(filter_dir, sample.id, ".", threshold, ".pass")
  hold <- paste0(filter_dir, sample.id, ".", threshold, ".hold")
  record <- paste0(filter_dir, sample.id, ".", threshold, ".record")
  jenks_estimate <- paste0(estimate_dir, sample.id, ".", threshold, ".jenks_estimate")

  # py4 <- paste0(python.script.dir, "/py4.filter.filter_with_pileup.py")
  # py4.command <- paste(py4, variant, intermediate, pass, hold, record,
  #                      jenks_estimate, threshold, SNP.database, depth)
  # system2(command = python, args = py4.command)

  filter_with_pileup_main(variant, intermediate, pass, hold, record, jenks_estimate,
                          threshold, SNP.database, depth)
  system2(command = "rm", args = variant)

}
