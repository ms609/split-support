# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/10_simulate.R")
}

on.exit(unlink("*.tmp.tre"))
for (aln in alnIDs) {
  if (file.exists(TNTFile(sim, aln, "ew"))) {
    message("Results found for ", aln)
  } else {
    system2(
      tntExec,
      paste0("run data-raw/tnt-ew.run", " ",
             DataFile(sim, aln), " ",
             TNTFile(sim, aln, "ew"), " ;")
    )
  }
}

# Validation step
for (aln in alnIDs) {
  thisFile <- TNTFile(sim, aln, "ew")
  if (file.exists(thisFile) && length(readLines(thisFile)) < 60) {
    message("Deleting incomplete analysis: ", aln)
    file.remove(thisFile)
  }
}
