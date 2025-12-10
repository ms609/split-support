# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/10_simulate.R")
}

session <- ssh::ssh_connect(hpcServer) # May prompt for authentication
remoteUser <- ssh::ssh_exec_internal(session, "echo $USER")$stdout |>
  rawToChar() |>
  trimws()
remoteWD <- file.path("", "nobackup", remoteUser, "spl-sup")
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
    
    remoteFile <- file.path(remoteWD, basename(MBFile(sim, aln)))
    writeLines(c(readLines("data-raw/slurm.sh"), paste("mb", remoteFile)),
               MBFile(aln.sh))
    
    ssh::ssh_exec_wait(session, paste("mkdir", "-p", remoteWD))
    ssh::scp_upload(session, MBFile(sim, aln), to = remoteWD)
    ssh::scp_upload(session, MBFile(aln.sh))
    
    ssh::ssh_exec_wait(session,
                       paste0("tr -d '\015' <", aln.sh, " >", aln.sh, ".sh"))
    print(ssh::ssh_exec_wait(session, paste0("sbatch ", aln.sh, ".sh")))
    
    # Remove unneeded results files
    outFiles <- list.files(path = mbDir, pattern = aln, full.names = TRUE)
    unlink(outFiles[-grep(paste0("(", paste0(keepExt, collapse = "|"), ")$"),
                          outFiles)])
  }
}

ssh::ssh_disconnect(session)
