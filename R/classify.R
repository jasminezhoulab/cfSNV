

extract_reads <- function(plasma.unmerged, reference,
                          samtools.dir, picard.dir, sample.id, python, java) {

  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  deeplearn_dir <- paste0(tmp.dir, "deeplearn/")
  if (system.file("extdata/tmp/deeplearn", package = "cfSNV") == "") {
    system2(command = "mkdir", args = deeplearn_dir)
  }

  ## DEDUPLICATE AND GENERATE BED FILE -----------------------------------------------
  filter_dir <- paste0(tmp.dir, "filter/")
  system2(command = "cat", args = paste0(
    filter_dir, sample.id, ".*.pass > ", deeplearn_dir, sample.id, ".pass"))
  system2(command = "cat", args = paste0(
    filter_dir, sample.id, ".*.hold > ", deeplearn_dir, sample.id, ".hold"))
  pass <- paste0(deeplearn_dir, sample.id, ".pass")
  hold <- paste0(deeplearn_dir, sample.id, ".hold")
  bed <- paste0(deeplearn_dir, sample.id, ".bed")

  py5 <- paste0(python.script.dir, "/py5.machinelearn.merge_dedup_candidates.py")
  py5.command <- paste(py5, pass, hold, bed, paste0(reference, ".fai"))
  system2(command = python, args = py5.command)

  system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.pass" ))
  system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.hold" ))
  system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.intermediate" ))
  system2(command = "rm", args = paste0(filter_dir, sample.id, ".*.record" ))
  system2(command = "rm", args = hold)
  system2(command = "rm", args = pass)

  ## EXTRACT READS AT BED FILE LOCATIONS ---------------------------------------------
  reference_dict <- paste0(deeplearn_dir, sample.id, '.reference.dict')
  interval_list <- paste0(deeplearn_dir, sample.id, ".interval_list")
  paired.reads.bam <- paste0(deeplearn_dir, sample.id, ".paired-reads.bam")
  paired.reads.qsort.bam <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.bam")
  paired.reads.qsort.sam <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.sam")
  log <- paste0(deeplearn_dir, sample.id, ".log")

  system2(command = java, args = paste(
    "-jar", picard.dir, "CreateSequenceDictionary", "-R", reference, "-O", reference_dict, ">>", log, "2>&1"))
  system2(command = java, args = paste0(
    "-Xms5G -jar ", picard.dir, " BedToIntervalList -I ", bed,
    " -O ", interval_list, " -SD ", reference_dict, " >> ", log, " 2>&1"))
  system2(command = java, args = paste0(
    "-Xms5G -jar ", picard.dir, " FilterSamReads -I ", plasma.unmerged,
    " -O ", paired.reads.bam, " -INTERVAL_LIST ", interval_list,
    " -FILTER includePairedIntervals -SORT_ORDER coordinate", " >> ", log, " 2>&1"))

  system2(command = samtools.dir, args = paste(
    "sort -n -o", paired.reads.qsort.bam, paired.reads.bam))
  system2(command = samtools.dir, args = paste(
    "view", paired.reads.qsort.bam, ">", paired.reads.qsort.sam))

  system2(command = "rm", args = reference_dict)
  system2(command = "rm", args = paired.reads.bam)
  system2(command = "rm", args = paired.reads.qsort.bam)
  system2(command = "rm", args = interval_list)
  system2(command = "rm", args = log)

}

prepare_bed_file <- function(reference, window=3, bedtools.dir, sample.id, python) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  deeplearn_dir <- paste0(tmp.dir, "deeplearn/")
  bed <- paste0(deeplearn_dir, sample.id, ".bed")
  getfastabed <- paste0(deeplearn_dir, sample.id, ".getfastabed")
  fastabed <- paste0(deeplearn_dir, sample.id, ".fastabed")
  preparebed <- paste0(deeplearn_dir, sample.id, ".preparebed")

  py6 <- paste0(python.script.dir, "/py6.machinelearn.bed_for_getfasta.py")
  py6.command <- paste(py6, window, bed, getfastabed)
  system2(command = python, args = py6.command)

  bedlog <- paste0(deeplearn_dir, sample.id, ".bedtools.log")
  system2(command = bedtools.dir, args = paste(
    "getfasta -bedOut -fi", reference, "-bed", getfastabed, ">", fastabed))

  system2(command = "rm", args = getfastabed)

  # py7 <- paste0(python.script.dir, "/py7.machinelearn.bed_for_feature.py")
  # py7.command <- paste(py7, window, fastabed, preparebed)
  # system2(command = python, args = py7.command)

}

extract_features_from_reads <- function(sample.id, python, window) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  deeplearn_dir <- paste0(tmp.dir, "deeplearn/")
  preparebed <- paste0(deeplearn_dir, sample.id, ".preparebed")
  paired.reads.qsort.sam <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.sam")
  paired.reads.qsort.features <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.features")
  filter_cluster.bed <- paste0(deeplearn_dir, sample.id, ".filter_cluster.bed")
  overlap.features <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.overlap.features")
  nonoverlap.features <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.nonoverlap.features")
  fastabed <- paste0(deeplearn_dir, sample.id, ".fastabed")

  # py8 <- paste0(python.script.dir, "/py8.machinelearn.extract_features_from_reads_filter_cluster.py")
  # py8.command <- paste(py8, preparebed, paired.reads.qsort.sam,
  #                      paired.reads.qsort.features, filter_cluster.bed)
  # system2(command = python, args = py8.command)
  machinelearn_extract_features_from_reads_filter_cluster(deeplearn_dir, deeplearn_dir, sample.id, window)

  system2(command = "rm", args = paired.reads.qsort.sam)
  system2(command = "rm", args = fastabed)

  py9 <- paste0(python.script.dir, "/py9.machinelearn.split_overlap_and_nonoverlap.py")
  py9.command <- paste(py9, paired.reads.qsort.features, overlap.features, nonoverlap.features)
  system2(command = python, args = py9.command)

  system2(command = "rm", args = paired.reads.qsort.features)

  py10.overlap <- paste0(python.script.dir, "/py10.machinelearn.expand_features_overlap.py")
  py10.overlap.command <- paste(py10.overlap, sample.id, overlap.features)
  system2(command = python, args = py10.overlap.command)

  py10.nonoverlap <- paste0(python.script.dir, "/py10.machinelearn.expand_features_nonoverlap.py")
  py10.nonoverlap.command <- paste(py10.nonoverlap, sample.id, nonoverlap.features)
  system2(command = python, args = py10.nonoverlap.command)

}

classify_and_generate_results <- function(bedtools.dir, sample.id, reference, python,
                                          MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  deeplearn_dir <- paste0(tmp.dir, "deeplearn/")
  overlap.features <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.overlap.features")
  nonoverlap.features <- paste0(deeplearn_dir, sample.id, ".paired-reads.qsort.nonoverlap.features")
  overlap.features.expand <- paste0(sample.id, ".paired-reads.qsort.overlap.features.expand")
  nonoverlap.features.expand <- paste0(sample.id, ".paired-reads.qsort.nonoverlap.features.expand")
  RF.overlap.model.gz <- paste0(extdata.dir, "/RF41_overlap.sklearn24.sav.gz")
  RF.nonoverlap.model.gz <- paste0(extdata.dir, "/RF41_nonoverlap.sklearn24.sav.gz")
  RF.overlap.model <- paste0(extdata.dir, "/RF41_overlap.sklearn24.sav")
  RF.nonoverlap.model <- paste0(extdata.dir, "/RF41_nonoverlap.sklearn24.sav")
  bed <- paste0(deeplearn_dir, sample.id, ".bed")

  system2(command = "gzip", args = paste("-cd", RF.overlap.model.gz, ">", RF.overlap.model))
  system2(command = "gzip", args = paste("-cd", RF.nonoverlap.model.gz, ">", RF.nonoverlap.model))

  ## CLASSIFY READS --------------------------------------------------------------------
  py11 <- paste0(python.script.dir, "/py11.machinelearn.RandomForest.py")
  py11.command <- paste(py11, deeplearn_dir, overlap.features.expand,
                        nonoverlap.features.expand, RF.overlap.model, RF.nonoverlap.model)
  system2(command = python, args = py11.command)

  ## GENERATE RESULTS ------------------------------------------------------------------
  py12 <- paste0(python.script.dir, "/py12.machinelearn.collect_results.py")
  py12.command <- paste(py12, deeplearn_dir, sample.id)
  system2(command = python, args = py12.command)

  system2(command = "rm", args = overlap.features)
  system2(command = "rm", args = nonoverlap.features)
  system2(command = "rm", args = paste0(deeplearn_dir, overlap.features.expand))
  system2(command = "rm", args = paste0(deeplearn_dir, nonoverlap.features.expand))
  system2(command = "rm", args = paste0(deeplearn_dir, overlap.features.expand, '.RFpred.csv'))
  system2(command = "rm", args = paste0(deeplearn_dir, nonoverlap.features.expand, '.RFpred.csv'))
  system2(command = "rm", args = bed)

  after_machine_learn <- paste0(deeplearn_dir, sample.id, ".after_machine_learn")
  after_indel <- paste0(deeplearn_dir, sample.id, ".after_indel")
  after_cluster <- paste0(deeplearn_dir, sample.id, ".after_cluster")

  rebase_dir <- paste0(tmp.dir, "rebase/")
  filter_indel.bed <- paste0(rebase_dir, sample.id, ".filter_indel.bed")
  filter_cluster.bed <- paste0(deeplearn_dir, sample.id, ".filter_cluster.bed")

  system2(command = bedtools.dir, args = paste(
    "intersect -v -a", after_machine_learn, "-b", filter_indel.bed, ">", after_indel))
  system2(command = bedtools.dir, args = paste(
    "intersect -v -a", after_indel, "-b", filter_cluster.bed, ">", after_cluster))

  system2(command = "rm", args = after_machine_learn)
  system2(command = "rm", args = after_indel)
  system2(command = "rm", args = filter_cluster.bed)
  system2(command = "rm", args = filter_indel.bed)

  py13 <- paste0(python.script.dir, "/py13.machinelearn.format_results.py")
  py13.command <- paste(py13, deeplearn_dir, sample.id, extdata.dir, paste0(reference, ".fai"),
                        MIN_HOLD_SUPPORT_COUNT, MIN_PASS_SUPPORT_COUNT)
  system2(command = python, args = py13.command)

  system2(command = "rm", args = after_cluster)
  system2(command = "rm", args = paste0(filter_cluster.bed, ".var.txt"))

  variant.report <- paste0(extdata.dir, "/", sample.id, ".variant_report.txt")
  jenks.estimate <- paste0(extdata.dir, "/", sample.id, ".jenks_estimate")

  variant.list <- read.delim(file = variant.report, header = FALSE)
  col_headings <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
  names(variant.list) <- col_headings

  estimate <- readLines(jenks.estimate, n = 1)
  tumor.fraction <- paste0(as.numeric(strsplit(estimate, '\t')[[1]][1])*100, "%")

  results <- list(variant.list = variant.list,
                  tumor.fraction = tumor.fraction)

  system2(command = "rm", args = variant.report)
  system2(command = "rm", args = jenks.estimate)
  return(results)
}


