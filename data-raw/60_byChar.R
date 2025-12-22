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
# 
# concordN <- vapply(cli::cli_progress_along(seq_len(nAln), "Analysing"),
#                    function (i) {
#   aln <- alnIDs[[i]]
#   
#   # Calculate concordances
#   dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(sim, aln))), 
#                                    nrow = nTip, byrow = TRUE,
#                                    dimnames = list(tips, NULL)))
#   
#   concCache <- ConcFile(sim, aln, "_chrN")
#   if (file.exists(concCache)) {
#     conc <- scan(concCache, quiet = TRUE)
#     if (length(conc) != nChar) {
#       file.remove(ConcFile(sim, aln))
#       stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
#     }
#   } else {
#     conc <- ClusteringConcordance(refSplits, dataset, normalize = TRUE,
#                                   return = "char")
#     write(conc, concCache)
#   }
#   conc
# }, double(nChar))

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
  rci = as.vector(consist[, "rci", ]),
  ri = as.vector(consist[, "ri", ]),
  rhi = 1 - as.vector(consist[, "rhi", ]),
  rhiBar = 1 - as.vector(consist[, "rhiBar", ]),
  conc = as.vector(concord),
  cat = rep(signif(cats, 3), each = nChar / nCats)
)

nona <- consDF |> na.omit()

CIPlot <- function(x) {
  boxplot(consDF[[x]] ~ consDF$cat, notch = TRUE,
          frame.plot = FALSE, las = 3,
          ylab = c(ci = "Consistency index",
                   ri = "Retention index",
                   rhi = "1 - Relative homoplasy index (median)",
                   rhiBar = "1 - Relative homoplasy index (mean)",
                   conc = "Clustering concordance",
                   rci = "Rescaled consistency index")[[x]],
          col = switch(x, conc = "gold", NULL),
          xlab = "Generative rate")
  cor.test(nona[[x]], nona$cat, method = "kendall")
}


par(mfrow = c(1, 5))

CIPlot("conc")    # tau ~ -0.6530646
CIPlot("ci")      # tau ~ -0.7734901
CIPlot("ri")      # tau ~ -0.4762535
CIPlot("rhi")     # tau ~ -0.5549066
#CIPlot("rhiBar")  # tau ~ -0.5431158
CIPlot("rci")     # tau ~ -0.5211873
