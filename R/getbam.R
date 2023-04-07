#' @title
#' Get Bam - align
#'
#' @description
#' Map raw reads from fastq (gzip or not) files to the reference genome.
#'
#' @export
getbam_align <- function(fastq1, fastq2, reference, SNP.database,
                    samtools.dir, picard.dir, bedtools.dir, GATK.dir, bwa.dir, sample.id,
                    output.dir, java.dir="java") {

  # time_label <- toString(Sys.time())
  # sample.id <- paste0(sample.id, "_", strsplit(time_label, " ")[[1]][2])

  ### Change the path of provided files
  # if (SNP.database == "dbSNP") {
  #   SNP.database = "/Users/huran/Desktop/cfSNV_development/demo/chr22_dbSNP.vcf"
  # }

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")
  if (system.file("extdata/tmp", package = "cfSNV") == "") {
    system2(command = "mkdir", args = tmp.dir)
  }
  reference_dict <- paste0(tmp.dir, sample.id, '.reference.dict')
  log <- paste0(tmp.dir, sample.id, ".log")
  sorted.bam <- paste0(tmp.dir, sample.id, ".sorted.bam")
  sorted.bam.bai <- paste0(tmp.dir, sample.id, ".sorted.bam.bai")
  sorted.rmdup.bam <- paste0(tmp.dir, sample.id, ".sorted.rmdup.bam")
  marked_dup_metrics.txt <- paste0(tmp.dir, sample.id, ".marked_dup_metrics.txt")
  sorted.rmdup.addreadgroup.bam <- paste0(tmp.dir, sample.id, ".sorted.rmdup.addreadgroup.bam")
  sorted.rmdup.addreadgroup.bam.bai <- paste0(tmp.dir, sample.id, ".sorted.rmdup.addreadgroup.bam.bai")
  forIndelRealigner.intervals <- paste0(tmp.dir, sample.id, ".forIndelRealigner.intervals")
  realignedBam.bam <- paste0(tmp.dir, sample.id, ".realignedBam.bam")
  realignedBam.bai <- paste0(tmp.dir, sample.id, ".realignedBam.bai")
  recal_data.table <- paste0(tmp.dir, sample.id, ".recal_data.table")
  recal.bam <- paste0(output.dir, "/", sample.id, ".recal.bam")
  intervals <- paste0(tmp.dir, sample.id, "*.intervals")
  table <- paste0(tmp.dir, sample.id, "*.table")
  txt <- paste0(tmp.dir, sample.id, "*.marked_dup_metrics.txt")

  ### REFERENCE DICT PREPARATION USING PICARD
  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " CreateSequenceDictionary REFERENCE=", reference, " OUTPUT=", reference_dict))

  system2(command = bwa.dir, args = paste("mem -t 4", reference, fastq1, fastq2, "|", samtools.dir,
                                      "sort -@3 -O BAM -o", sorted.bam, "-"))

  system2(command = samtools.dir, args = paste("index", sorted.bam))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " MarkDuplicates I=", sorted.bam, " O=", sorted.rmdup.bam,
    " M=", marked_dup_metrics.txt, " REMOVE_DUPLICATES=true"))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " AddOrReplaceReadGroups I=", sorted.rmdup.bam, " O=", sorted.rmdup.addreadgroup.bam,
    " RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"))

  system2(command = samtools.dir, args = paste("index", sorted.rmdup.addreadgroup.bam))

  system2(command = "rm", args = paste(sorted.bam, sorted.bam.bai, sorted.rmdup.bam))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T RealignerTargetCreator -R", reference, "-I", sorted.rmdup.addreadgroup.bam,
    "-o", forIndelRealigner.intervals))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T IndelRealigner -R", reference, "-I", sorted.rmdup.addreadgroup.bam,
    "-targetIntervals", forIndelRealigner.intervals, "-o", realignedBam.bam))

  system2(command = "rm", args = paste(sorted.rmdup.addreadgroup.bam.bai, sorted.rmdup.addreadgroup.bam))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T BaseRecalibrator -R", reference, "-I", realignedBam.bam,
    "-o", recal_data.table, "-knownSites", SNP.database))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T PrintReads -R", reference, "-I", realignedBam.bam,
    "-BQSR", recal_data.table, "-o", recal.bam))

  system2(command = "rm", args = paste(reference_dict))
  system2(command = "rm", args = paste(realignedBam.bam, realignedBam.bai))
  system2(command = "rm", args = paste(intervals, table, txt))

}

#' @title
#' Get Bam - align after merge
#'
#' @description
#' First merge the overlapping read mates in cfDNA sequencing data and
#' then map raw reads from fastq (gzip or not) files to the reference genome
#'
#' @export
getbam_align_after_merge <- function(fastq1, fastq2, reference, SNP.database,
                                     samtools.dir, picard.dir, bedtools.dir,
                                     GATK.dir, bwa.dir, flash.dir, sample.id,
                                     output.dir, java.dir="java") {

  # time_label <- toString(Sys.time())
  # sample.id <- paste0(sample.id, "_", strsplit(time_label, " ")[[1]][2])

  ### Change the path of provided files
  # if (SNP.database == "dbSNP") {
  #   SNP.database = "/Users/huran/Desktop/cfSNV_development/demo/chr22_dbSNP.vcf"
  # }

  extdata.dir <- system.file("extdata", package = "cfSNV", mustWork = TRUE)
  tmp.dir <- paste0(extdata.dir, "/tmp/")
  if (system.file("extdata/tmp", package = "cfSNV") == "") {
    system2(command = "mkdir", args = tmp.dir)
  }
  reference_dict <- paste0(tmp.dir, sample.id, '.reference.dict')
  log <- paste0(tmp.dir, sample.id, ".log")
  notCombined_1.fastq.gz <- paste0(tmp.dir, sample.id, '.notCombined_1.fastq.gz')
  notCombined_2.fastq.gz <- paste0(tmp.dir, sample.id, '.notCombined_2.fastq.gz')
  notCombined.sorted.bam <- paste0(tmp.dir, sample.id, '.notCombined.sorted.bam')
  notCombined.sorted.bam.bai <- paste0(tmp.dir, sample.id, '.notCombined.sorted.bam.bai')
  extendedFrags.fastq.gz <- paste0(tmp.dir, sample.id, '.extendedFrags.fastq.gz')
  extendedFrags.sorted.bam <- paste0(tmp.dir, sample.id, '.extendedFrags.sorted.bam')
  extendedFrags.sorted.bam.bai <- paste0(tmp.dir, sample.id, '.extendedFrags.sorted.bam.bai')
  extendedFrags.sorted.rmdup.bam <- paste0(tmp.dir, sample.id, '.extendedFrags.sorted.rmdup.bam')
  extendedFrags.marked_dup_metrics.txt <- paste0(tmp.dir, sample.id, '.extendedFrags.marked_dup_metrics.txt')
  extendedFrags.sorted.rmdup.addreadgroup.bam <- paste0(tmp.dir, sample.id, '.extendedFrags.sorted.rmdup.addreadgroup.bam')
  extendedFrags.sorted.rmdup.addreadgroup.bam.bai <- paste0(tmp.dir, sample.id, '.extendedFrags.sorted.rmdup.addreadgroup.bam.bai')
  notCombined.sorted.rmdup.bam <- paste0(tmp.dir, sample.id, '.notCombined.sorted.rmdup.bam')
  notCombined.marked_dup_metrics.txt <- paste0(tmp.dir, sample.id, '.notCombined.marked_dup_metrics.txt')
  notCombined.sorted.rmdup.addreadgroup.bam <- paste0(tmp.dir, sample.id, '.notCombined.sorted.rmdup.addreadgroup.bam')
  notCombined.sorted.rmdup.addreadgroup.bam.bai <- paste0(tmp.dir, sample.id, '.notCombined.sorted.rmdup.addreadgroup.bam.bai')
  extendedFrags.forIndelRealigner.intervals <- paste0(tmp.dir, sample.id, '.extendedFrags.forIndelRealigner.intervals')
  extendedFrags.realignedBam.bam <- paste0(tmp.dir, sample.id, '.extendedFrags.realignedBam.bam')
  extendedFrags.realignedBam.bai <- paste0(tmp.dir, sample.id, '.extendedFrags.realignedBam.bai')
  notCombined.forIndelRealigner.intervals <- paste0(tmp.dir, sample.id, '.notCombined.forIndelRealigner.intervals')
  notCombined.realignedBam.bam <- paste0(tmp.dir, sample.id, '.notCombined.realignedBam.bam')
  notCombined.realignedBam.bai <- paste0(tmp.dir, sample.id, '.notCombined.realignedBam.bai')
  extendedFrags.recal_data.table <- paste0(tmp.dir, sample.id, '.extendedFrags.recal_data.table')
  extendedFrags.recal.bam <- paste0(output.dir, "/", sample.id, '.extendedFrags.recal.bam')
  notCombined.recal_data.table <- paste0(tmp.dir, sample.id, '.notCombined.recal_data.table')
  notCombined.recal.bam <- paste0(output.dir, "/", sample.id, '.notCombined.recal.bam')

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " CreateSequenceDictionary REFERENCE=", reference, " OUTPUT=", reference_dict))

  system2(command = flash.dir, args = paste("-M 20 -x 0.01 -z -q -d", tmp.dir, "-o", sample.id, fastq1, fastq2))

  system2(command = bwa.dir, args = paste("mem -t 4", reference, notCombined_1.fastq.gz, notCombined_2.fastq.gz,
                                      "|", samtools.dir, "sort -@3 -O BAM -o", notCombined.sorted.bam, "-"))

  system2(command = bwa.dir, args = paste("mem -t 4", reference, extendedFrags.fastq.gz,
                                      "|", samtools.dir, "sort -@3 -O BAM -o", extendedFrags.sorted.bam, "-"))

  system2(command = samtools.dir, args = paste("index", extendedFrags.sorted.bam))

  system2(command = samtools.dir, args = paste("index", notCombined.sorted.bam))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " MarkDuplicates I=", extendedFrags.sorted.bam,
    " O=", extendedFrags.sorted.rmdup.bam, " M=", extendedFrags.marked_dup_metrics.txt, " REMOVE_DUPLICATES=true"))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " AddOrReplaceReadGroups I=", extendedFrags.sorted.rmdup.bam,
    " O=", extendedFrags.sorted.rmdup.addreadgroup.bam, " RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"))

  system2(command = "rm", args = paste(extendedFrags.sorted.rmdup.bam, extendedFrags.sorted.bam,
                                       extendedFrags.sorted.bam.bai))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " MarkDuplicates I=", notCombined.sorted.bam,
    " O=", notCombined.sorted.rmdup.bam, " M=", notCombined.marked_dup_metrics.txt, " REMOVE_DUPLICATES=true"))

  system2(command = java.dir, args = paste0(
    "-Xms1g -jar ", picard.dir, " AddOrReplaceReadGroups I=", notCombined.sorted.rmdup.bam,
    " O=", notCombined.sorted.rmdup.addreadgroup.bam, " RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"))

  system2(command = samtools.dir, args = paste("index", notCombined.sorted.rmdup.addreadgroup.bam))

  system2(command = samtools.dir, args = paste("index", extendedFrags.sorted.rmdup.addreadgroup.bam))

  system2(command = "rm", args = paste(notCombined.sorted.rmdup.bam, notCombined.sorted.bam,
                                       notCombined.sorted.bam.bai))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T RealignerTargetCreator", "-R", reference, "-I", extendedFrags.sorted.rmdup.addreadgroup.bam,
    "-o", extendedFrags.forIndelRealigner.intervals))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T IndelRealigner", "-R", reference, "-I", extendedFrags.sorted.rmdup.addreadgroup.bam,
    "-targetIntervals", extendedFrags.forIndelRealigner.intervals, "-o", extendedFrags.realignedBam.bam))

  system2(command = "rm", args = paste(extendedFrags.sorted.rmdup.addreadgroup.bam,
                                       extendedFrags.sorted.rmdup.addreadgroup.bam.bai))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T RealignerTargetCreator", "-R", reference, "-I", notCombined.sorted.rmdup.addreadgroup.bam,
    "-o", notCombined.forIndelRealigner.intervals))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T IndelRealigner", "-R", reference, "-I", notCombined.sorted.rmdup.addreadgroup.bam,
    "-targetIntervals", notCombined.forIndelRealigner.intervals, "-o", notCombined.realignedBam.bam))

  system2(command = "rm", args = paste(notCombined.sorted.rmdup.addreadgroup.bam,
                                       notCombined.sorted.rmdup.addreadgroup.bam.bai))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T BaseRecalibrator", "-R", reference, "-I", extendedFrags.realignedBam.bam,
    "-o", extendedFrags.recal_data.table, "-knownSites", SNP.database))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T PrintReads", "-R", reference, "-I", extendedFrags.realignedBam.bam,
    "-BQSR", extendedFrags.recal_data.table, "-o", extendedFrags.recal.bam))

  system2(command = "rm", args = paste(extendedFrags.realignedBam.bam,
                                       extendedFrags.realignedBam.bai))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T BaseRecalibrator", "-R", reference, "-I", notCombined.realignedBam.bam,
    "-o", notCombined.recal_data.table, "-knownSites", SNP.database))

  system2(command = java.dir, args = paste(
    "-Xms1g -jar", GATK.dir, "-T PrintReads", "-R", reference, "-I", notCombined.realignedBam.bam,
    "-BQSR", notCombined.recal_data.table, "-o", notCombined.recal.bam))

  system2(command = "rm", args = paste(notCombined.realignedBam.bam,
                                       notCombined.realignedBam.bai))

  intervals <- paste0(tmp.dir, sample.id, "*.intervals")
  table <- paste0(tmp.dir, sample.id, "*.table")
  txt <- paste0(tmp.dir, sample.id, "*.marked_dup_metrics.txt")

  hist <- paste0(tmp.dir, sample.id, ".hist*")

  system2(command = "rm", args = paste(intervals, table, txt))

  system2(command = "rm", args = paste(reference_dict, extendedFrags.fastq.gz,
                                       notCombined_1.fastq.gz, notCombined_2.fastq.gz, hist))
}
