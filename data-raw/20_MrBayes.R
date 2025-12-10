source("data-raw/config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/simulate.R")
}


template <- readLines("data-raw/mb.nex")

for (aln in alns) {
  
  if (file.exists(MBFile(aln, "tstat"))) {
    message("Tree probabilities found for alignment ", aln)
  } else {
    on.exit(unlink(MBFile(aln)))
    writeLines(
      c(readLines(DataFile(aln)), template),
      MBFile(aln)
    )
    system2(mbExec, MBFile(aln))
    
    
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
