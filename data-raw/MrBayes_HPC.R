source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

ssh::ssh_connect(hpcServer) # Will prompt for authentication

template <- readLines("split-support/mb.nex")

for (aln in alns) {
  
  if (file.exists(MBFile(aln, "tstat"))) {
    message("Tree probabilities found for alignment ", aln)
  } else {
    aln.sh <- paste0(aln, ".sh")
    on.exit(unlink(c(MBFile(aln), MBFile(aln.sh))))
    writeLines(
      c(readLines(DataFile(aln)), template),
      MBFile(aln)
    )
    
    writeLines(c(readLines("split-support/slurm.sh"),
                 paste("mb", MBFile(aln))),
               MBFile(aln.sh))
    
    scp_upload(session, MBFile(aln))
    scp_upload(session, MBFile(aln.sh))
    
    print(ssh_exec_wait(session, paste0("tr -d '\015' <", aln.sh, " >", aln.sh, ".sh")))
    print(ssh_exec_wait(session, paste0("sbatch ", aln.sh, ".sh")))
    
    # Download the file back and verify it is the same
    scp_download(session, "slurm-1536761.out", to = tempdir())
    if (tools::md5sum(file_path) != 
        tools::md5sum(file.path(tempdir(), "slurm-1536761.out"))) {
      warning("MD5 mismatch: ", aln)
    }
    
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

ssh::ssh_disconnect(session)