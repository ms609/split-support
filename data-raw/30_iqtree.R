# Set simulation identifier
sim <- "gam"

source("data-raw/_config.R")

if(!dir.exists("data-raw/alignments")) {
  source("data-raw/10_simulate.R")
}

for (aln in alnIDs) {
  if (file.exists(IQFile(sim, aln, ".splits.nex"))) {
    message("Results found for ", aln)
  } else {
    seqs <- ape::read.nexus.data(DataFile(sim, aln))
    phyle <- IQFile(sim, aln)
    on.exit(unlink(phyle))
    writeLines(c(
      paste(length(seqs), length(seqs[[1]])),
      "",
      apply(cbind(names(seqs),
            vapply(seqs, paste0, character(1), collapse = "")), 1, paste,
            collapse = "\t")
      ), phyle)
    system2(
      iqExec,
      # http://www.iqtree.org/doc/Command-Reference
      paste0(" -s ", phyle,
             " -st DNA ", # Sequence type: DNA
                        # Model: Jukes-Cantor / + 6 gamma rate cats
             " -mset ", switch(sim, "aln" = "JC", "gam" = "JC+G6",
                               stop("unknown `sim`")),
             "  -mrate E ", # Equal rates only
             " -nt 1", # Number of threads; -nt auto -ntmax 6 is slow
             " -seed 1 ", # Set random seed for reproducibility
             #" -b 1000", # Nonparametric bootstrap is slow and can't
                          # be run alongside UFB
             " -bb 1000 ", # Number of ultrafast bootstrap replicates
             " -bnni ", # Avoids branch support overestimates in UFB
             " -lbp 1000 ", # Fast local bootstrap probability
             " -alrt 1000 ", # Approximate likelihood ratio test
             " -abayes ", # Approximate Bayes test
             " -quiet ",
             " --redo-tree ", # Overwrite previous run results
             ""
             )
    )
    
    # Remove unneeded results files
    iqKeepExt <- c(
      "treefile", # Maximum likelihood tree
      #"iqtree", # Lists order of support values
      # Here: SH-aLRT support (%) / local bootstrap support (%) / 
      #       aBayes support / ultrafast bootstrap support (%)
      #"contree", # UF Bootstrap consensus tree
      "splits\\.nex" # UF Bootstrap split supports
    )
    
    outFiles <- list.files(path = iqDir, pattern = aln, full.names = TRUE)
    
    unlink(outFiles[-grep(paste0("(", paste0(iqKeepExt, collapse = "|"), ")$"),
                          outFiles)])
  }
}
