### Analytical setup ###

# Size of tree and dataset
nTip <- 48
tips <- paste0("t", seq_len(nTip))
nChar <- nTip * 2
nCats <- 6

# Number of alignments to analyse
nAln <- 1000
alnIDs <- formatC(seq_len(nAln), width = 4, flag = 0)

### Location of executables ###

# Download command-driven TNT from 
# http://www.lillo.org.ar/phylogeny/tnt/ZIPCHTNT.ZIP
# and unzip to a convenient local folder.

# Command to launch TNT executable
tntExec <- "C:/Programs/Phylogeny/tnt/tnt.exe"

# Download MrBayes from 
# https://nbisweden.github.io/MrBayes/download.html
# and unzip/install to a convenient local folder.

# Command to launch MrBayes executable
mbExec <- "C:/Programs/Phylogeny/MrBayes/bin/mb.3.2.7-win64.exe"

# Download IQ-tree from http://www.iqtree.org/
# and unzip/install to a convenient local folder.

# Command to launch IQ-tree executable
iqExec <- "C:/Programs/Phylogeny/iqtree-3.0.1-Windows/bin/iqtree3.exe"

# Username and server ID for HPC login
# Ideally, set a system environment variable, perhaps by running
# file.edit("~/.Renviron")
# and adding a line
# sshLogin = username@hpc.institution.ac.uk
# Then restart the R environment
hpcServer <- Sys.getenv("sshLogin")
# Alternatively, set directly with 
# hpcServer <- "username@hpc.institution.ac.uk"

### Location of output files ###
iqDir <- "data-raw/iqtree/"
mbDir <- "data-raw/MrBayes/"
tntDir <- "data-raw/TNT/"
alnDir <- "data-raw/alignments/"
concDir <- "data-raw/concordance/"
hDir <- "data-raw/entropy/"

# Set up directory structure
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Patterns to use when creating files
CreateDir(concDir)
ConcFile <- function(sim, id, suffix = "") {
  file.path(concDir, paste0(sim, id, suffix, ".txt"))
}

CreateDir(hDir)
EntropyFile <- function(sim, id) {
  file.path(hDir, paste0(sim, id, ".txt"))
}

PartQFile <- function(sim, id) {
  file.path(hDir, paste0(sim, id, ".q.txt"))
}

CreateDir(alnDir)
DataFile <- function(sim, id, ext = ".nex") {
  file.path(alnDir, paste0(sim, id, ext))
}

CreateDir(iqDir)
IQFile <- function(sim, id, suffix = "") {
  file.path(iqDir, paste0(sim, id, ".phy", suffix))
}

CreateDir(mbDir)
MBFile <- function(sim, id = "", suffix = NULL) {
  file.path(mbDir, paste0(sim, id, if(!is.null(suffix)) ".", suffix))
}

CreateDir(tntDir)
TNTFile <- function(sim, id, wt = "ew") {
  file.path(tntDir, paste0(sim, id, ".", wt, ".out"))
}

# MrBayes output files to retain
keepExt <- c(
  "con\\.tre", # Consensus tree - why not
  "parts", "tstat", # Partitions and probabilities
  # "trprobs" # Sampled trees and probabilities
  # "mcmc" # Standard deviations of splits - see `tstat`
  "pstat" # Convergence diagnostics
)

Panel <- function(i, xOffset = 3, yOffset = 0) {
  if (is.numeric(i)) {
    i <- paste0("(", letters[i], ")")
  }
  usr <- par("usr")
  xOffset <- xOffset * strwidth("M")
  yOffset <- yOffset * strheight("M")
  
  x <- usr[[1]] - xOffset
  if (par("xlog")) {
    x <- 10 ^ x
  }
  
  y <- usr[[4]] - yOffset
  if (par("ylog")) {
    y <- 10 ^ y
  }
  
  text(x, y, i, xpd = NA, adj = c(1, 1)) # adj: right-align | top-align
}
