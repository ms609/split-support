# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/10_simulate.R")
}


template <- readLines(file.path("data-raw", sprintf("mb-%s.nex", sim)))

for (aln in alnIDs) {
  
  if (file.exists(MBFile(sim, aln, "tstat"))) {
    message("Tree probabilities found for alignment ", aln)
  } else {
    on.exit(unlink(MBFile(sim, aln)))
    writeLines(
      c(readLines(DataFile(sim, aln)), template),
      MBFile(sim, aln)
    )
    system2(mbExec, MBFile(sim, aln))
    
    
    # Remove unneeded results files
    keepExt <- c(
      "con\\.tre", # Consensus tree - why not
      "parts", "tstat", # Partitions and probabilities
      # "trprobs" # Sampled trees and probabilities
      # "mcmc" # Standard deviations of splits - see `tstat`
      "pstat" # Convergence diagnostics
    )
    
    outFiles <- list.files(path = mbDir, pattern = aln, full.names = TRUE)
    
    unlink(outFiles[-grep(paste0("(", paste0(keepExt, collapse = "|"), ")$"),
                          outFiles)])
    
  }
}
