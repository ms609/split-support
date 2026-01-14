# Load required libraries
library("TreeTools")
library("TreeSearch")

# Load configuration settings
source("data-raw/_config.R")

# Set simulation identifier here
sim <- "gam"
cats <- phangorn::discrete.gamma(1, nCats)

referenceTree <- file.path("data-raw", sprintf("reference-%s.tre", sim)) |>
  read.tree()
refSplits <- as.Splits(referenceTree)
tips <- names(read.nexus.data(DataFile(sim, "0001")))

concord <- vapply(
  cli::cli_progress_along(seq_len(nAln), "Analysing"),
  function(i) {
    aln <- alnIDs[[i]]
    
    # Calculate concordances
    dataset <- DataFile(sim, aln) |>
      read.nexus.data() |>
      unlist() |>
      matrix(nrow = nTip, byrow = TRUE,
             dimnames = list(tips, NULL)) |>
      MatrixToPhyDat()
    
    concCache <- ConcFile(sim, aln, "_chr")
    if (file.exists(concCache)) {
      conc <- scan(concCache, quiet = TRUE)
      if (length(conc) != nChar) {
        file.remove(concCache)
        stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
      }
    } else {
      conc <- ClusteringConcordance(refSplits, dataset, normalize = FALSE,
                                    return = "char")
      write(conc, concCache)
    }
    conc
  },
  double(nChar)
)

qConcord <- vapply(
  cli::cli_progress_along(seq_len(nAln), "Analysing"),
  function(i) {
    aln <- alnIDs[[i]]
    
    # Calculate concordances
    dataset <- DataFile(sim, aln) |>
      read.nexus.data() |>
      unlist() |>
      matrix(nrow = nTip, byrow = TRUE,
             dimnames = list(tips, NULL)) |>
      MatrixToPhyDat()
    
    concCache <- ConcFile(sim, aln, "_chrQ")
    if (file.exists(concCache)) {
      conc <- scan(concCache, quiet = TRUE)
      if (length(conc) != nChar) {
        file.remove(concCache)
        stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
      }
    } else {
      conc <- QuartetConcordance(refSplits, dataset, return = "char")
      write(conc, concCache)
    }
    conc
  },
  double(nChar)
)

consist <- vapply(
  cli::cli_progress_along(seq_len(nAln), "Analysing"),
  function(i) {
    aln <- alnIDs[[i]]
    
    # Calculate concordances
    dataset <- DataFile(sim, aln) |>
      read.nexus.data() |>
      unlist() |>
      matrix(nrow = nTip, byrow = TRUE,
             dimnames = list(tips, NULL)) |>
      MatrixToPhyDat()
    
    consCache <- ConcFile(sim, aln, "_cns")
    if (file.exists(consCache)) {
      cons <- do.call(cbind, read.table(consCache))
      if (nrow(cons) != nChar) {
        file.remove(ConcFile(sim, aln, "_cns"))
        stop("Dimension mismatch; is consistency cache ", aln, " out of date?")
      }
    } else {
      cons <- TreeSearch::Consistency(dataset, referenceTree, nRelabel = 1000)
      write.table(cons, consCache)
    }
    cons
  },
  matrix(double(), nChar, 7)
)

# Compile statistics into single data.frame for analysis and filtering
consDF <- data.frame(
  ci = as.vector(consist[, "ci", ]),
  rci = as.vector(consist[, "rci", ]),
  ri = as.vector(consist[, "ri", ]),
  rhi = 1 - as.vector(consist[, "rhi", ]),
  rhiBar = 1 - as.vector(consist[, "rhiBar", ]),
  conc = as.vector(concord),
  qConc = as.vector(qConcord),
  cat = rep(sprintf("%.3f", cats), each = nChar / nCats)
)

nona <- consDF |> na.omit()

# Calculate resolution:
resolution <- apply(nona, 2, function(x) length(unique(x)))
resolution[rev(order(resolution)[-1])]

# Plot relationship with generative rate categories
CIPlot <- function(x, calcTau = FALSE) {
  boxplot(
    consDF[[x]] ~ consDF$cat,
    frame.plot = FALSE,
    las = 3,
    range = 1,
    pch = 16,
    cex = 0.2,
    ylab = c(
      ci = "Consistency index",
      ri = "Retention index",
      rhi = "1 - Relative homoplasy index",
      rhiBar = "1 - Relative homoplasy index (mean)",
      conc = "Clustering concordance",
      qConc = "Quartet concordance",
      rci = "Rescaled consistency index"
    )[[x]],
    col = switch(x, conc = "gold", NULL),
    xlab = ""
  )
  
  if (calcTau) {
    cor.test(nona[[x]], as.numeric(nona$cat), method = "kendall")
  }
}

{
  pdf("Fig 3 - character concordance.pdf", 8.4, 2.4)
  par(
    mfrow = c(1, 6),
    mar = c(3.7, 5.2, 0.2, 0.2),
    cex = 0.6,
    oma = c(1, 0, 0, 0)
  )

  CIPlot("conc", calcTau = TRUE) # tau ~ -0.6530646
  CIPlot("qConc", calcTau = TRUE) # tau ~ -0.4571022
  CIPlot("ci", calcTau = TRUE) # tau ~ -0.7734901
  CIPlot("ri", calcTau = TRUE) # tau ~ -0.4762535
  CIPlot("rhi", calcTau = TRUE) # tau ~ -0.5549066
  CIPlot("rci", calcTau = TRUE) # tau ~ -0.5211873
  mtext("Generative rate for character", 1, line = 0, outer = TRUE, cex = 0.6)
  dev.off()
}
