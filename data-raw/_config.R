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
iqExec <- "C:/Programs/Phylogeny/iqtree-2.2.0-Windows/bin/iqtree2.exe"

hpcServer <- "pjjg18@hamilton8.dur.ac.uk"

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
ConcFile <- function(sim, id) {
  file.path(concDir, paste0(sim, id, ".txt"))
}

# Patterns to use when creating files
CreateDir(hDir)
EntropyFile <- function(sim, id) {
  file.path(hDir, paste0(sim, id, ".txt"))
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
