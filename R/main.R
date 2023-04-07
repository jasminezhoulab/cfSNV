#' @title
#' Parameter Recommendation
#'
#' @description
#' Report the range of sequencing depth, the roughly estimated tumor fraction, and the detection limit range.
#'
#' @param plasma.unmerged a BAM file of the cfDNA sequencing data.
#' @param normal a BAM file of normal sample sequencing data.
#' @param plasma.merged.extendedFrags a BAM file of extendedFrags.
#' @param plasma.merge.notCombined a BAM file of reads not combined.
#' @param target.bed a BED file of target regions.
#' @param reference a FASTA file of reference genome.
#' @param SNP.database a VCF file of SNP database.
#' @param samtools.dir the path to an executable Samtools.
#' @param sample.id a sample name.
#' @param roughly_estimated_tf if calculate the roughly estimated tumor fraction.
#' @param python.dir the path to Python.
#'
#' @examples
#' ## input files
#' demo.dir <- '/Users/huran/Desktop/cfSNV_development/demo/'
#' plasma.unmerged <- paste0(demo.dir, 'plasma.recal.bam')
#' normal <- paste0(demo.dir, 'normal.recal.bam')
#' plasma.merged.extendedFrags <- paste0(demo.dir, 'plasma.extendedFrags.recal.bam')
#' plasma.merge.notCombined <- paste0(demo.dir, 'plasma.notCombined.recal.bam')
#' target.bed <- paste0(demo.dir, 'example_target_regions.bed')
#' reference <- paste0(demo.dir, 'chr22.fa')
#' SNP.database <- paste0(demo.dir, 'chr22_dbSNP.vcf')
#' samtools.dir <- '/usr/local/bin/samtools'
#' sample.id <- '1st'
#' python.dir <- '/usr/local/bin/python3'
#'
#'
#' parameter_recommend(
#'   plasma.unmerged, normal,
#'   plasma.merged.extendedFrags, plasma.merge.notCombined,
#'   target.bed, reference, SNP.database, samtools.dir, sample.id, python.dir=python.dir
#' )
#' @export
parameter_recommend <- function(plasma.unmerged, normal,
                            plasma.merged.extendedFrags, plasma.merge.notCombined,
                            target.bed, reference, SNP.database, samtools.dir,
                            sample.id, roughly_estimated_tf = FALSE, python.dir) {

  # time_label <- toString(Sys.time())
  # sample.id <- paste0(sample.id, "_", strsplit(time_label, " ")[[1]][2])

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")
  estimate_dir <- paste0(tmp.dir, "estimate/")
  if (system.file("extdata/tmp", package = "cfSNV") == "") {
    system2(command = "mkdir", args = tmp.dir)
  }

  MIN_PASS_SUPPORT_COUNT <- 5
  SNP.database.processed <- paste0(tmp.dir, sample.id, ".dbSNP.vcf")

  ### Change the path of provided files
  # if (SNP.database == "dbSNP") {
  #   SNP.database = "/Users/huran/Desktop/cfSNV_development/demo/chr22_dbSNP.vcf"
  #   system2(command = "grep", args = paste("grep SNV", SNP.database, "| grep -Eiv \"GENEINFO\" >", SNP.database.processed))
  # } else if (SNP.database == "None") {
  #   SNP.database.processed = "/Users/huran/Desktop/cfSNV_development/demo/empty.vcf"
  # } else { # user provide file path
  #   system2(command = "grep", args = paste("-Eiv \"#\"", SNP.database, ">", SNP.database.processed))
  # }

  system2(command = "grep", args = paste("-Eiv \"#\"", SNP.database, ">", SNP.database.processed))

  SNP.database <- SNP.database.processed

  depth_dir <- paste0(tmp.dir, "depth/")
  if (system.file("extdata/tmp/depth", package = "cfSNV") == "") {
    system2(command = "mkdir", args = depth_dir)
  }

  # system2(command = "wget", args = paste0("--load-cookies /tmp/cookies.txt ",
  #                                        "\"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies ",
  #                                        "--no-check-certificate 'https://docs.google.com/uc?export=download&id=1FwdXIr4keuKROA_gNX0OSrOcgLajVn7n' ",
  #                                        "-O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=1FwdXIr4keuKROA_gNX0OSrOcgLajVn7n\" -O ",
  #                                        tmp.dir, "dbSNP.vcf && rm -rf /tmp/cookies.txt")
  # )

  if (roughly_estimated_tf) {
    depth <- get_sequencing_depth(plasma.unmerged, target.bed, samtools.dir, sample.id, python.dir)

    ## data preparation
    pileup_and_rebase(plasma.unmerged, normal, plasma.merged.extendedFrags,
                      plasma.merge.notCombined, target.bed, reference,
                      samtools.dir, depth, sample.id, python.dir, MIN_PASS_SUPPORT_COUNT)

    collect_indel_related_regions(sample.id, reference, python.dir)

    ## iteration
    threshold <- 1.0000

    while (TRUE) {
      iterate_filter(SNP.database, depth, threshold, sample.id, python.dir, MIN_PASS_SUPPORT_COUNT)

      jenks_estimate <- paste0(estimate_dir, sample.id, ".", threshold, ".jenks_estimate")
      estimate_results <- readLines(jenks_estimate, n = 1)

      fields <- strsplit(estimate_results, '\t')[[1]]
      n_pass_hold <- as.double(fields[1])
      threshold <- as.double(fields[3])
      estimate <- as.double(fields[2])

      system2(command = "rm", args = jenks_estimate)

      if (n_pass_hold == 0 | threshold == 0) {
        break
      }
      break
    }
    rebase_dir <- paste0(tmp.dir, "rebase/")
    filter_dir <- paste0(tmp.dir, "filter/")
    indel <- paste0(rebase_dir, sample.id, ".indel")
    rebase <- paste0(rebase_dir, sample.id, ".rebase")
    filter_indel.bed <- paste0(rebase_dir, sample.id, ".filter_indel.bed")
    system2(command = "rm", args = indel)
    system2(command = "rm", args = rebase)
    system2(command = "rm", args = filter_indel.bed)
    system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.pass" ))
    system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.hold" ))
    system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.intermediate" ))
    system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.record" ))
  }

  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)
  # reticulate::source_python(paste0(python.script.dir, "/py15.estimate_sequencing_depth.py"))

  depth_dir <- paste0(tmp.dir, "depth/")
  region_depth <- paste0(depth_dir, sample.id, ".region_depth.txt")
  region_depth_output <- paste0(depth_dir, sample.id, ".region_depth_out.txt")
  system2(command = samtools.dir,
          args = paste0("bedcov ", target.bed, " ", plasma.unmerged,
                        " | awk \'BEGIN {FS=OFS=", "\"", "\t", "\"","}",
                        " {print $4/($3-$2)}\' > ", region_depth))

  # depth_results <- depth_range(region_depth)
  py16 <- paste0(python.script.dir, "/py16.estimate_depth_range.py")
  py16.command <- paste(py16, region_depth, region_depth_output)
  system2(command = python.dir, args = py16.command)

  depth_file <- file(region_depth_output, "r")
  avrg <- as.numeric(readLines(depth_file, n = 1))
  med <- as.numeric(readLines(depth_file, n = 1))
  p95 <- as.numeric(readLines(depth_file, n = 1))
  close(depth_file)

  cat("The per base coverage of the plasma sample for each genomic region in the target bed file:\n")
  cat(paste0("average = ", signif(avrg, 6)))
  cat(paste0(", median = ", signif(med, 6)))
  # p95 <- as.numeric(depth_results[3])
  cat(paste0(", 95th percentile = ", signif(p95, 6)), '\n')

  if (roughly_estimated_tf) {
    cat(paste0("\nThe roughly estimated tumor fraction in the plasma sample: ", signif(estimate*100, 6), "%"), "\n")
    cat("For a more accurate estimation, please run DetectMuts.", "\n")
  }

  if ((signif(5/p95, 4)<=1) & (signif(12/p95, 4)<=1)) {
    cat(paste0("\nLowest detectable VAF range under the default parameters: [", signif(5/p95, 4)*100, "%, ", signif(12/p95, 4)*100, "%]"),"\n")

    cat("\nTo detect different levels of lowest VAF,", "\n")
    if (p95 >= 500) {
      cat(paste0("at 0.5% VAF: minHold = ",
                 round(0.005*p95, digits = 0)+6, ", minPass = ",
                 round(0.005*p95, digits = 0), ";"), "\n")
    }
    cat(paste0("at 1% VAF: minHold = ",
               round(0.01*p95, digits = 0)+6, ", minPass = ",
               round(0.01*p95, digits = 0), ";"), "\n")
    cat(paste0("at 5% VAF: minHold = ",
               round(0.05*p95, digits = 0)+6, ", minPass = ",
               round(0.05*p95, digits = 0)), "\n")
    cat("Note: decreasing the parameters (i.e. minHold and minPass) \ncan lower the detection limit, but may also lower the variant quality.")
  } else {
    cat("The per base coverage is too low for high-quality mutation calling.", "\n")
  }

  system2(command = "rm", args = region_depth_output)
  system2(command = "rm", args = SNP.database.processed)
  system2(command = "rm", args = region_depth)
  #system2(command = "rm", args = paste("-r", tmp.dir))
}



#' @title
#' Variant Calling
#'
#' @description
#' This function allows you to call somatic mutations.
#'
#' @param plasma.unmerged a BAM file of the cfDNA sequencing data.
#' @param normal a BAM file of normal sample sequencing data.
#' @param plasma.merged.extendedFrags a BAM file of extendedFrags.
#' @param plasma.merge.notCombined a BAM file of reads not combined.
#' @param target.bed a BED file of target regions.
#' @param reference a .fa file of reference genome.
#' @param SNP.database a .vcf of SNP database.
#' @param samtools.dir the path to an executable Samtools.
#' @param picard.dir the path to picard.jar.
#' @param bedtools.dir the path to an executable BEDTools.
#' @param sample.id the name of sample.
#' @param MIN_HOLD_SUPPORT_COUNT default is 12.
#' @param MIN_PASS_SUPPORT_COUNT default is 5.
#' @param java.dir the path to java, default is "java".
#' @param python.dir the path to Python.
#'
#' @return a variant list and a tumor fraction
#'
#' @examples
#' ## input files
#' demo.dir <- '/Users/huran/Desktop/cfSNV_development/demo/'
#' plasma.unmerged <- paste0(demo.dir, 'plasma.recal.bam')
#' normal <- paste0(demo.dir, 'normal.recal.bam')
#' plasma.merged.extendedFrags <- paste0(demo.dir, 'plasma.extendedFrags.recal.bam')
#' plasma.merge.notCombined <- paste0(demo.dir, 'plasma.notCombined.recal.bam')
#' target.bed <- paste0(demo.dir, 'example_target_regions.bed')
#' reference <- paste0(demo.dir, 'chr22.fa')
#' SNP.database <- paste0(demo.dir, 'chr22_dbSNP.vcf')
#' samtools.dir <- '/usr/local/bin/samtools'
#' picard.dir <- '/usr/local/bin/picard.jar'
#' bedtools.dir <- '/usr/local/bin/bedtools2/bin/bedtools'
#' sample.id <- '1st'
#' MIN_HOLD_SUPPORT_COUNT <- 9
#' MIN_PASS_SUPPORT_COUNT <- 3
#' python.dir <- '/usr/local/bin/python3'
#'
#' results <- variant_calling(
#'   plasma.unmerged, normal, plasma.merged.extendedFrags, plasma.merge.notCombined,
#'   target.bed, reference, SNP.database, samtools.dir, picard.dir, bedtools.dir,
#'   sample.id, MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT, python.dir=python.dir
#' )
#' @export
variant_calling <- function(plasma.unmerged, normal,
                            plasma.merged.extendedFrags, plasma.merge.notCombined,
                            target.bed, reference, SNP.database,
                            samtools.dir, picard.dir, bedtools.dir, sample.id,
                            MIN_HOLD_SUPPORT_COUNT = 12, MIN_PASS_SUPPORT_COUNT = 5,
                            java.dir = "java", python.dir) {

  # time_label <- toString(Sys.time())
  # sample.id <- paste0(sample.id, "_", strsplit(time_label, " ")[[1]][2])

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")
  estimate_dir <- paste0(tmp.dir, "estimate/")
  if (system.file("extdata/tmp", package = "cfSNV") == "") {
    system2(command = "mkdir", args = tmp.dir)
  }

  SNP.database.processed <- paste0(tmp.dir, sample.id, ".dbSNP.vcf")

  ### Change the path of provided files
  # if (SNP.database == "dbSNP") {
  #   SNP.database = "/Users/huran/Desktop/cfSNV_development/demo/chr22_dbSNP.vcf"
  #   system2(command = "grep", args = paste("grep SNV", SNP.database, "| grep -Eiv \"GENEINFO\" >", SNP.database.processed))
  # } else if (SNP.database == "None") {
  #   SNP.database.processed = "/Users/huran/Desktop/cfSNV_development/demo/empty.vcf"
  # } else { # user provide file path
  #   system2(command = "grep", args = paste("-Eiv \"#\"", SNP.database, ">", SNP.database.processed))
  # }

  system2(command = "grep", args = paste("-Eiv \"#\"", SNP.database, ">", SNP.database.processed))

  SNP.database <- SNP.database.processed

  depth <- get_sequencing_depth(plasma.unmerged, target.bed, samtools.dir, sample.id, python.dir)

  ## data preparation
  pileup_and_rebase(plasma.unmerged, normal, plasma.merged.extendedFrags,
                    plasma.merge.notCombined, target.bed, reference,
                    samtools.dir, depth, sample.id, python.dir, MIN_PASS_SUPPORT_COUNT)

  collect_indel_related_regions(sample.id, reference, python.dir)

  ## iteration
  threshold <- 1.0000

  while (TRUE) {
    iterate_filter(SNP.database, depth, threshold, sample.id, python.dir, MIN_PASS_SUPPORT_COUNT)

    jenks_estimate <- paste0(estimate_dir, sample.id, ".", threshold, ".jenks_estimate")
    estimate_results <- readLines(jenks_estimate, n = 1)

    fields <- strsplit(estimate_results, '\t')[[1]]
    n_pass_hold <- as.double(fields[1])
    threshold <- as.double(fields[3])

    system2(command = "rm", args = jenks_estimate)

    if (n_pass_hold == 0 | threshold == 0) {
      break
    }
  }

  rebase_dir <- paste0(tmp.dir, "rebase/")
  indel <- paste0(rebase_dir, sample.id, ".indel")
  rebase <- paste0(rebase_dir, sample.id, ".rebase")
  system2(command = "rm", args = indel)
  system2(command = "rm", args = rebase)

  ## classification
  extract_reads(plasma.unmerged, reference, samtools.dir, picard.dir, sample.id, python.dir, java.dir)

  window <- 3
  prepare_bed_file(reference, window, bedtools.dir, sample.id, python.dir)

  extract_features_from_reads(sample.id, python.dir, window)

  results <- classify_and_generate_results(bedtools.dir, sample.id, reference, python.dir,
                                MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT)

  system2(command = "rm", args = SNP.database.processed)
  #system2(command = "rm", args = paste("-r", tmp.dir))
  return(results)
}
