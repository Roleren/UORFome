#' Create STAR genome index
#'
#' Used as reference when aligning data
STAR.index <- function(GTF, fasta, output.dir,
                       rrna = NULL, trna = NULL, ncrna = NULL, phix = NULL,
                       script = "/export/valenfs/projects/Pipelines/STAR_Aligner/STAR_MAKE_INDEX.sh") {
  phix <- ifelse(is.null(phix), "", paste("-phix", phix))
  rrna <- ifelse(is.null(rrna), "", paste("-rRna", rrna))
  trna <- ifelse(is.null(trna), "", paste("-tRNA", trna))
  ncrna <- ifelse(is.null(ncrna), "", paste("-ncRNA", ncrna))
  GTF <- paste("-genomeGTF", GTF)
  fasta <- paste("-genome", fasta)

  full <- paste(script, "-o", output.dir, phix, rrna, trna, ncrna, GTF, fasta)
  print(full)
  if (.Platform$OS.type == "unix") {
    system(command = full)
    message("Alignment done")
  } else stop("STAR is not supported on windows!")
}

#' Align data with STAR
STAR.folder <- function(input.dir, output.dir, index.dir,
                        paired.end = "no",
                        steps = "tr-ge", adapter.sequence = "auto",
                        min.length = 15, trim.front = 3,
                        alignment.type = "Local", max.cpus = 90,
                        include.subfolders = "n",
                        script = "/export/valenfs/projects/Pipelines/STAR_Aligner/RNA_Align_pipeline_folder.sh -h") {

  full <- paste(script, "-f", input.dir, "-o", output.dir, "-p", paired.end,
                "-l", min.length, "-g", index.dir, "-s", steps,
                "-a", adapter.sequence, "-t", trim.front,
                "-A", alignment.type, "-m", max.cpus, "-S", include.subfolders)
  print(full)
  if (.Platform$OS.type == "unix") {
    system(command = full)
    message("Alignment done")
  } else stop("STAR is not supported on windows!")
}

STAR.merge.tcp <- function(input.dir) {
  dir.create(paste0(input.dir, "/merged"))
  system("samtools --help")
}
