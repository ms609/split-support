# Load required libraries
library("TreeTools")
library("TreeSearch")

# Load configuration settings
source("data-raw/_config.R")

# Set simulation identifier here
sim <- "gam"

referenceTree <- file.path("data-raw", sprintf("reference-%s.tre", sim)) |>
  read.tree()
refSplits <- as.Splits(referenceTree)
tips <- names(read.nexus.data(DataFile(sim, "0001")))

concordN <- numeric(0)
for (i in cli::cli_progress_along(seq_len(nAln), "Analysing")) {
  aln <- alnIDs[[i]]
  
  # Calculate concordances
  dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(sim, aln))), 
                                   nrow = nTip, byrow = TRUE,
                                   dimnames = list(tips, NULL)))
  
  concCache <- ConcFile(sim, aln, "_chrN")
  if (file.exists(concCache)) {
    conc <- scan(concCache, quiet = TRUE)
    if (length(conc) != nChar) {
      file.remove(ConcFile(sim, aln))
      stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
    }
  } else {
    conc <- ClusteringConcordance(refSplits, dataset, normalize = TRUE,
                                  return = "char")
    write(conc, concCache)
  }
  concordN <- cbind(concordN, conc)
}

concN <- data.frame(
  conc = as.vector(concordN),
  cat = rep(signif(cats, 4), each = nChar / nCats)
)
boxplot(concN$conc ~ concN$cat, notch = TRUE,
        frame.plot = FALSE, las = 3,
        ylab = "Normalized clustering concordance",
        xlab = "Generative rate")
#plot(conc~log(cat), data = concN)
summary(lm(conc ~ cat, data = concN))
summary(lm(conc ~ log(cat), data = concN))
