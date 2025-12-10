# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/10_simulate.R")
}

session <- ssh::ssh_connect(hpcServer) # May prompt for authentication
template <- readLines(file.path("data-raw", sprintf("mb-%s.nex", sim)))

for (aln in alnIDs) {
  
  if (file.exists(MBFile(sim, aln, "tstat"))) {
    message("Tree probabilities found for alignment ", aln)
  } else {
    aln.sh <- paste0(sim, aln, ".sh")
    on.exit(unlink(c(MBFile(sim, aln), MBFile(aln.sh))))
    writeLines(
      c(readLines(DataFile(sim, aln)), template),
      MBFile(sim, aln)
    )
    
    writeLines(c(readLines("data-raw/slurm.sh"),
                 paste("mb", MBFile(sim, aln))),
               MBFile(aln.sh))
    
    ssh::scp_upload(session, MBFile(sim, aln))
    ssh::scp_upload(session, MBFile(aln.sh))
    
    print(ssh::ssh_exec_wait(session, paste0("tr -d '\015' <", aln.sh, " >", aln.sh, ".sh")))
    print(ssh::ssh_exec_wait(session, paste0("sbatch ", aln.sh, ".sh")))
    
    # Download the file back and verify it is the same
    ssh::scp_download(session, "slurm-1536761.out", to = tempdir())
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
