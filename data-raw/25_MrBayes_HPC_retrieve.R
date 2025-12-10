# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

session <- ssh::ssh_connect(hpcServer) # May prompt for authentication
remoteUser <- ssh::ssh_exec_internal(session, "echo $USER")$stdout |>
  rawToChar() |>
  trimws()
remoteWD <- file.path("", "nobackup", remoteUser, "spl-sup")

for (aln in alnIDs) {
  
  if (file.exists(MBFile(sim, aln, "tstat"))) {
    message("Tree probabilities found locally for alignment ", aln)
  } else {
    remoteFiles <- gsub("\\", "", fixed = TRUE,
                        paste0(file.path(remoteWD, basename(MBFile(sim, aln))),
                               ".", keepExt))
    
    if (all(vapply(remoteFiles, function(f) {
      ssh::ssh_exec_wait(session, sprintf("test -e %s", shQuote(f))) == 0
    }, logical(1)))) {
      for (rf in remoteFiles) {
        status <- ssh::scp_download(session, rf, to = MBFile(""))
      }
    } else {
      message("No results for ", sim, aln)
    }
  }
}

ssh::ssh_disconnect(session)
