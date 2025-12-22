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

concordN <- vapply(cli::cli_progress_along(seq_len(nAln), "Analysing"),
                   function (i) {
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
  conc
}, double(nChar))

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
        file.remove(ConcFile(sim, aln, "_chr"))
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
        # file.remove(ConcFile(sim, aln, "_cns"))
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

consDF <- data.frame(
  ci = as.vector(consist[, "ci", ]),
  cci = as.vector(consist[, "cci", ]),
  ri = as.vector(consist[, "ri", ]),
  rhi = 1 - as.vector(consist[, "rhi", ]),
  rhiBar = 1 - as.vector(consist[, "rhiBar", ]),
  rci = as.vector(consist[, "rci", ]),
  cat = rep(signif(cats, 3), each = nChar / nCats)
)

CIPlot <- function(x) {
  boxplot(consDF[[x]] ~ consDF$cat, notch = TRUE,
          frame.plot = FALSE, las = 3,
          ylab = c(ci = "Consistency index",
                   cci = "CCI",
                   ri = "Retention index",
                   rhi = "1 - Relative homoplasy index (median)",
                   rhiBar = "1 - Relative homoplasy index (mean)",
                   rci = "Rescaled consistency index")[[x]],
          xlab = "Generative rate")
  cor.test(consDF[[x]], consDF$cat, method = "kendall") # -0.6690298 
}


par(mfrow = c(4, 2))

concDF <- data.frame(
  conc = as.vector(concord),
  cat = rep(signif(cats, 4), each = nChar / nCats)
)

boxplot(concDF$conc ~ concDF$cat, notch = TRUE,
        frame.plot = FALSE, las = 3,
        ylab = "Clustering concordance",
        xlab = "Generative rate")
cor.test(concDF$conc, concDF$cat, method = "kendall") # -0.6690298 


concN <- data.frame(
  conc = as.vector(concordN),
  cat = rep(signif(cats, 4), each = nChar / nCats)
)
boxplot(concN$conc ~ concN$cat, notch = TRUE,
        frame.plot = FALSE, las = 3,
        ylab = "Normalized clustering concordance",
        xlab = "Generative rate")
cor.test(concN$conc, concN$cat, method = "kendall") # tau = -0.5654185 

CIPlot("ci")      # tau ~ -0.7885806
CIPlot("cci")     # tau ~ -0.4323527
CIPlot("ri")      # tau ~ -0.4762535
CIPlot("rhi")     # tau ~ -0.5549066
CIPlot("rhiBar")  # tau ~ -0.5431158
CIPlot("rci")     # tau ~ -0.5211873
