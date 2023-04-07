

# utils::globalVariables(c("estimate_sequencing_depth", "_cfSNV_timesTwo"))

## Get Depth from BAM ---------------------------------------------------

get_sequencing_depth <- function(plasma.recal.bam, target.bed, samtools.dir, sample.id, python) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")
  if (system.file("extdata/tmp", package = "cfSNV") == "") {
    system2(command = "mkdir", args = tmp.dir)
  }

  depth_dir <- paste0(tmp.dir, "depth/")
  if (system.file("extdata/tmp/depth", package = "cfSNV") == "") {
    system2(command = "mkdir", args = depth_dir)
  }

  base.mapped <- paste0(depth_dir, sample.id, ".base_mapped.txt")

  ## call Samtools
  system2(command = samtools.dir, args = paste(
    "stats", plasma.recal.bam, "| grep -i 'bases mapped (cigar):' | cut -f3 >>", base.mapped))


  ## call py15
  output <- paste0(depth_dir, sample.id, ".depth.txt")
  # reticulate::source_python(paste0(python.script.dir, "/py15.estimate_sequencing_depth.py"))
  # depth <- estimate_sequencing_depth(target.bed, base.mapped, output)

  py15 <- paste0(python.script.dir, "/py15.estimate_sequencing_depth.py")
  py15.command <- paste(py15, target.bed, base.mapped, output)
  system2(command = python, args = py15.command)

  depth_file <- file(output, "r")
  depth <- as.numeric(readLines(depth_file, n = 1))
  close(depth_file)

  system2(command = "rm", args = base.mapped)
  system2(command = "rm", args = output)
  return(depth)
}

## Pileup and Rebase ---------------------------------------------------

pileup_and_rebase <- function(plasma.unmerged, normal,
                              plasma.merged.extendedFrags, plasma.merge.notCombined,
                              target.bed, reference, samtools.dir, depth, sample.id, python,
                              MIN_PASS_SUPPORT_COUNT) {
  python.script.dir <- system.file("python", package = "cfSNV", mustWork = TRUE)

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")

  pileup_dir <- paste0(tmp.dir, "pileup/")
  if (system.file("extdata/tmp/pileup", package = "cfSNV") == "") {
    system2(command = "mkdir", args = pileup_dir)
  }

  rebase_dir <- paste0(tmp.dir, "rebase/")
  if (system.file("extdata/tmp/rebase", package = "cfSNV") == "") {
    system2(command = "mkdir", args = rebase_dir)
  }

  ## call Samtools mpileup
  pileup.file <- paste0(pileup_dir, sample.id, ".pileup")
  log <- paste0(pileup_dir, sample.id, ".log")

  system2(command = samtools.dir, args = paste(
    "mpileup -A -q 0 -Q 3 -B --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -s -f",
    reference, "-l", target.bed, plasma.unmerged, normal,
    plasma.merged.extendedFrags, plasma.merge.notCombined, ">", pileup.file, "2>", log))

  ## call py1 to get rebase file
  output <- paste0(rebase_dir, sample.id, ".rebase")
  indel <- paste0(rebase_dir, sample.id, ".indel")

  # py1 <- paste0(python.script.dir, "/py1.dataprep.rebase.quicker.py")
  # py1.command <- paste(py1, pileup.file, output, indel, depth)
  # system2(command = python, args = py1.command)

  ParseFile(pileup.file, output, indel, depth, MIN_PASS_SUPPORT_COUNT)
  system2(command = "rm", args = pileup.file)
  system2(command = "rm", args = log)

}
