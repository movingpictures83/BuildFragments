#context("alpine")

source("plugins/BuildFragments/core.R")
source("plugins/BuildFragments/fit_bias.R")
source("plugins/BuildFragments/helper.R")
source("plugins/BuildFragments/vlmm.R")
source("plugins/BuildFragments/plots.R")
source("plugins/BuildFragments/estimate_abundance.R")
source("plugins/BuildFragments/predict.R")
library(alpineData)
  library(GenomicAlignments)
  library(rtracklayer)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.NCBI.GRCh38)


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  pfix <<- prefix()
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
}
run <- function() {}

output <- function(outputfile) {
  gap <- ERR188088()
  bam.file <- c("ERR188088" = file.path(paste(pfix, parameters["bamfile", 2], sep="/")))
  export(gap, con=bam.file)
  load(file=paste(pfix, parameters["preproc", 2], sep="/"))
  readlength <- as.integer(parameters["readlength", 2])
  minsize <- as.integer(parameters["minsize", 2])
  maxsize <- as.integer(parameters["maxsize", 2])
  gene.names <- names(ebt.fit)[6:8]
  names(gene.names) <- gene.names
  fragtypes <- lapply(gene.names, function(gene.name) {
    buildFragtypes(ebt.fit[[gene.name]],
                   Hsapiens, readlength,
                   minsize, maxsize)
  })
  saveRDS(fragtypes, outputfile)
}
