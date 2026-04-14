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

partCorrectList <- vector("list", nAln)
partQualList <- vector("list", nAln)
postProbList <- vector("list", nAln)
concordList <- vector("list", nAln)
bremerList <- vector("list", nAln)
splitHList <- vector("list", nAln)
tntStatList <- vector("list", nAln)
iqStatList <- vector("list", nAln)

tntStats <- c("symFq", "symGC", "boot", "jak", "pois")
iqStats <- c("alrt", "lbp", "abayes", "ufb", "sCF") # .iqtree output file gives order

for (i in cli::cli_progress_along(seq_len(nAln), "Analysing")) {
  aln <- alnIDs[[i]]
  
  # Load MrBayes partitions
  parts <- read.table(MBFile(sim, aln, "parts"), skip = 2 + nTip)
  partitions <- setNames(as.Splits(parts[, 2], tips), paste0("mb", parts[, 1]))
  
  # Load TNT partitions
  tntFile <- TNTFile(sim, aln, "ew")
  if (!file.exists(tntFile)) {
    warning("No TNT results available for ", aln, ".")
    next;
  }
  tntTree <- suppressMessages(ReadTntTree(tntFile, tipLabels = tips))
  if (!inherits(tntTree, "multiPhylo")) {
    warning("Only one tree found in file ", aln, "; missing TNT output.")
    next;
  }
  if (!all.equal(tntTree[[1]], tntTree[[2]])) {
    # Trees may differ in resolution: partitions with 0 Bremer support will
    # be collapsed in tntTree[[1]].
    if (any(!as.Splits(tntTree[[1]]) %in% as.Splits(tntTree[[2]]))) {
      warning("Trees don't match in ", aln, ". Check TNT output.")
      next;
    }
  }
  tntParts <- as.Splits(tntTree[[2]])
  nTntNode <- length(tntParts)
  tntOnly <- !tntParts %in% partitions
  if (any(tntOnly)) {
    partitions <- c(
      partitions,
      setNames(tntParts[[tntOnly]], paste0("tnt", seq_len(sum(tntOnly))))
    )
  }
  
  # Load IQ-tree partitions
  iqTreeFile <- IQFile(sim, aln, ".treefile")
  scfFile <- paste0(iqTreeFile, ".cf.tree")
  if (!file.exists(scfFile)) {
    df <- DataFile(sim, aln)
    dfLines <- readLines(df)
    dfLines[[5]] <- "FORMAT DATATYPE=DNA;"
    tmpFile <- paste0(df, ".tmp")
    file.rename(df, tmpFile)
    on.exit({
      if (file.exists(tmpFile)) {
        unlink(df)
        file.rename(tmpFile, df)
      }
    }, add = TRUE)
    writeLines(dfLines, df)
    system2(Sys.getenv("iqtree.exe"),
            c("-t", iqTreeFile, "-s", df, "--scf 100000", "-nt 4"),
            stdout = NULL)
    unlink(df)
    file.rename(tmpFile, df)
  }

  iqTree <- read.tree(scfFile)
  iqParts <- as.Splits(iqTree)
  iqOnly <- !iqParts %in% partitions
  if (any(iqOnly)) {
    partitions <- c(
      partitions,
      setNames(iqParts[[iqOnly]], paste0("iq", seq_len(sum(iqOnly))))
    )
  }
  ufbLines <- readLines(IQFile(sim, aln, ".splits.nex"))[-seq_len(nTip + 13)]
  ufbLines <- ufbLines[seq_len(which.max(ufbLines == ";") - 1)]
  ufbLines <- do.call(rbind, strsplit(trimws(ufbLines), "\t"))
  ufbParts <- as.Splits(t(vapply(
    strsplit(trimws(gsub(",", "", fixed = TRUE, ufbLines[, 2])), " "),
    function(a) {
      tabulate(as.numeric(a), nTip) == 1
    }, logical(nTip))), tipLabels = tips)
  trivial <- TrivialSplits(ufbParts)
  ufbVals <- as.numeric(ufbLines[!trivial, 1])
  ufbParts <- ufbParts[[!trivial]]
  ufbOnly <- !ufbParts %in% partitions
  if (any(ufbOnly)) {
    partitions <- c(
      partitions,
      setNames(ufbParts[[ufbOnly]], paste0("ufb", seq_len(sum(ufbOnly))))
    )
  }
  
  # Once all partitions are loaded, label where possible
  
  # Populate TNT supports
  brem <- rep(NA_real_, length(partitions))
  partId2 <- match(as.Splits(tntTree[[2]]), partitions)
  brem[partId2] <- 0
  partBrem <- tntTree[[1]]$node.label[
    as.numeric(names(as.Splits(tntTree[[1]]))) - NTip(tntTree[[1]])]
  brem[match(as.Splits(tntTree[[1]]), partitions)] <- as.numeric(partBrem)
  
  tags <- strsplit(tntTree[[2]]$node.label, "/")
  partTags <- tags[as.numeric(names(as.Splits(tntTree[[2]]))) - NTip(tntTree[[2]])]
  tntTags <- matrix(NA_real_, length(partitions), length(tntStats),
                    dimnames = list(NULL, tntStats))
  tntTags[partId2, ] <- t(vapply(partTags, function(tag) {
    x <- gsub("[", "-", fixed = TRUE, gsub("]", "", fixed = TRUE, tag))
    x[x == "?"] <- NA_real_
    x[x == "???"] <- Inf
    as.numeric(x)
  }, numeric(length(tntStats))))
  
  
  # Populate IQ-tree supports
  iqMatch <- match(iqParts, partitions)
  tags <- strsplit(iqTree[["node.label"]], "/")
  partTags <- tags[as.numeric(names(iqParts)) - NTip(iqTree)]
  iqTags <- matrix(NA_real_, length(partitions), length(iqStats),
                   dimnames = list(NULL, iqStats))
  iqTags[iqMatch, ] <- t(vapply(partTags, as.numeric, numeric(length(iqStats))))
  
  
  # Populate Ultra-Fast bootstrap supports for partitions not in consensus
  ufbMatch <- match(ufbParts, partitions)
  iqTags[ufbMatch, "ufb"] <- ufbVals
  
  
  # Populate posterior probabilities
  pp <- read.table(MBFile(sim, aln, "tstat"), skip = 1,
                   header = TRUE, comment.char = "")
  pp <- setNames(pp[, "Probability..s."], pp[, "ID"])
  
  
  # Calculate concordances
  dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(sim, aln))), 
                                   nrow = nTip, byrow = TRUE,
                                   dimnames = list(tips, NULL)))
  
  concCache <- ConcFile(sim, aln)
  if (file.exists(concCache)) {
    conc <- as.matrix(read.table(concCache))
    if (dim(conc)[[1]] != dim(tntTags)[[1]]) {
      file.remove(ConcFile(sim, aln))
      stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
    }
    if (!"wQuartet" %in% colnames(conc)) {
      conc <- cbind(quartet = QuartetConcordance(partitions, dataset,
                                                 weight = FALSE),
                    wQuartet = conc[, 1],
                    
                    conc[, -1]
                    )
      write.table(conc, concCache)
    }
  } else {
    # For efficiency, calculate the complete concordance statistics once and
    # derive the associated measures below.
    # Output matches that produced via ClusteringConcordance(normalize = T/F)
    cAll <- ClusteringConcordance(
      partitions,
      dataset,
      normalize = FALSE,
      return = "all"
    )
    bestSums <- rowSums(cAll["hBest", , ])
    .Rezero <- function(value, zero) {
      (value - zero) / (1 - zero)
    }

    conc <- cbind(
      quartet = QuartetConcordance(partitions, dataset, weight = FALSE),
      wQuartet = QuartetConcordance(partitions, dataset),
      cluster = rowSums(cAll["mi", , ]) / bestSums, # = ClustConc(norm = FALSE)
      phylo = PhylogeneticConcordance(partitions, dataset),
      mutual = MutualClusteringConcordance(partitions, dataset),
      # Not used in this study:
      shared = SharedPhylogeneticConcordance(partitions, dataset),
      clusterNorm = .Rezero(
        rowSums(cAll["mi", , ]) / bestSums,
        rowSums(cAll["miRand", , ]) / bestSums
      ) # = ClustConc(partitions, dataset, norm = TRUE)
    )
    write.table(conc, concCache)
  }
  
  hCache <- EntropyFile(sim, aln)
  if (file.exists(hCache)) {
    h <- as.matrix(read.table(hCache))
    if (dim(h)[[1]] != dim(tntTags)[[1]]) {
      file.remove(hCache)
      stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
    }
  } else {
    h <- cbind(
      clustering = TreeDist::ClusteringEntropy(partitions, sum = FALSE),
      # Not used in this study:
      splitwise = TreeDist::SplitwiseInfo(partitions, sum = FALSE)
    )
    write.table(h, hCache)
  }
  
  partInRef <- partitions %in% refSplits
  partCorrectList[[i]] <- partInRef
  
  qCache <- PartQFile(sim, aln)
  if (file.exists(qCache)) {
    partQ <- scan(qCache, quiet = TRUE)
  } else {
    partQ <- rep(1, length(partInRef))
    partQ[!partInRef] <- vapply(
      seq_along(partitions)[!partInRef],
      function(i) TreeDist::MutualClusteringInfo(partitions[[i]], refSplits,
                                                 normalize = h[, "clustering"][[i]]),
      double(1))
    write(partQ, qCache)
  }
  partQualList[[i]] <- partQ
  
  postProbList[[i]] <- c(pp, rep(0, sum(tntOnly, iqOnly, ufbOnly)))
  concordList[[i]] <- conc
  splitHList[[i]] <- h
  bremerList[[i]] <- brem
  tntStatList[[i]] <- tntTags
  iqStatList[[i]] <- iqTags
  
}; cli::cli_progress_done()

# Reformat lists into vectors/matrices
partCorrect <- do.call(c, partCorrectList)
partQual <- do.call(c, partQualList)
postProb <- do.call(c, postProbList)
concord <- do.call(rbind, concordList)
splitH <- do.call(rbind, splitHList)
bremer <- do.call(c, bremerList)
tntStat <- do.call(rbind, tntStatList)
iqStat <- do.call(rbind, iqStatList)

colnames(tntStat) <- tntStats
colnames(iqStat) <- iqStats

common <- rowSums(is.na(concord)) == 0 &
  rowSums(is.na(tntStat)) == 0 &
  !is.na(bremer) &
  rowSums(is.na(iqStat)) == 0

# Arrange in data.frame to allow subsequent filtering and analysis
allDat <- data.frame(
  occurs = partCorrect,
  partQual,
  postProb,
  tntStat,
  iqStat,
  bremer,
  concord
)

dat <- data.frame(occurs = partCorrect, partQual, postProb, concord) |> na.omit()

# Compute Somers' D, from which the C-index may be derived
SomersD <- function(score, target) {
  # C index = (Dxy + 1) / 2
  fit <- Hmisc::rcorr.cens(score, target)
  
  est <- fit["Dxy"]
  # Standard error for Dxy
  se  <- fit["S.D."] / sqrt(fit["n"]) 
  
  ci  <- est + c(-1, 1) * 1.96 * se
  list(estimate = est, ci95 = ci)
}
# Derive the C-index from Somers' D
CIndex <- function(score, target) {
  lapply(SomersD(score, target), function(x) (x + 1) / 2)
}


# How well does a measure predict whether a split is in the true tree?
# We set `cf` to include only splits for which data is available under
# both `var` and `cf`, to allow a straight comparison.
Histy <- function(var, breaks = 16, even = TRUE, cf = var) { # "Mosaic plot"
  entries <- !is.na(var) & !is.na(cf)
  outcomes <- factor(partCorrect[entries], levels = c("FALSE", "TRUE"),
                     ordered = TRUE)
  var <- var[entries]
  brks <- if (isTRUE(even)) {
    quantile(var, seq(0, 1, length.out = breaks))
  } else if (even == "log") {
    quantile(var, log(1:breaks) / log(breaks))
  } else {
    seq(min(var), max(var), length.out = breaks)
  }
  brks <- unique(brks)
  pattern <- if (max(brks) > 2) "%.0f" else "%.3f"
  binLabels <- sprintf(pattern, brks[-length(brks)])
  bins <- cut(var, breaks = brks)
  
  col <- c("TRUE" = "3", "FALSE" = "2")
  call <- match.call()
  if (substr(as.character(call[-1])[[1]], 1, 7) != "concord") {
    col <- adjustcolor(col, alpha.f = 0.5)
  }
  title <- as.character(call[2])
  title <- switch(
    title,
    "postProb" = "Posterior probability",
    "concord[, \"cluster\"]" = "Clustering concordance",
    "concord[, \"quartet\"]" = "Quartet concordance",
    "concord[, \"wQuartet\"]" = "Weighted quartet conc.",
    "concord[, \"mutual\"]" = "Mutual clustering concordance",
    "bremer" = "Bremer support",
    "tntStat[, \"symFq\"]" = "Symmetric frequency",
    "tntStat[, \"symGC\"]" = "Groups pres / cont",
    "tntStat[, \"boot\"]" = "TNT bootstrap",
    "tntStat[, \"jak\"]" = "TNT jackknife",
    "tntStat[, \"pois\"]" = "Poisson resampling",
    
    "iqStat[, \"ufb\"]" = "Ultra-fast bootstrap",
    "iqStat[, \"lbp\"]" = "Local bootstrap probabilities",
    "iqStat[, \"alrt\"]" = "Approx. lik. ratio test",
    "iqStat[, \"abayes\"]" = "Approx. Bayes",
    "iqStat[, \"sCF\"]" = "Site concordance factor",
    title
  )
  
  tab <- table(bins, outcomes)
  spTab <- spineplot(
    tab,
    main = title,
    col = col,
    axes = FALSE,
    xaxlabels = "",
    yaxlabels = "",
    xlab = "",
    ylab = "",
    border = NA
  )

  binCounts <- rowSums(tab)
  #binLabels <- binLabels[binCounts > 0]
  widths <- binCounts / sum(binCounts)
  edges <- cumsum(widths)
  leftEdges <- c(0, edges[-length(edges)])
  # The centre is the average of the left and right edges
  centres <- (leftEdges + edges) / 2
  usr <- par("usr")
  plotWidth <- usr[[2]] - usr[[1]]
  # Map our 0-1 centres to the actual usr coordinates
  x <- usr[[1]] + (centres * plotWidth)

  text(
    x = x,
    y = -0.06,
    labels = binLabels,
    srt = 90,
    adj = 1,
    xpd = NA,
    cex = 0.8
  )

  cacheFile <- file.path("data-raw", "roc", 
                         gsub("[ ,\"\\[]|\\]", "",
                              paste(call[-1], collapse = "-")))
  message(cacheFile)
  if (file.exists(cacheFile)) {
    load(cacheFile)
  } else {
    roc <- pROC::roc(predictor = var, response = as.numeric(outcomes),
                     quiet = TRUE)
    cIdx <- CIndex(var, partQual[entries])
    save(roc, cIdx, file = cacheFile)
  }

  message("n = ", sum(entries), ": ", title)

  mtext(
    paste0(
      "ROC-AUC = ",
      sprintf("%.2f", roc$auc),
      "; ",
      "C-index = ",
      sprintf("%.2f", cIdx$estimate)
    ),
    3,
    cex = 0.6
  )
}

# Produce figure as PDF
{
  cairo_pdf("Fig 2 - edge concordance.pdf", 5.4, 8.4)
  par(mar = c(1.6, 0.8, 3, 0.8), font.main = 1, cex.main = 0.9)
  yAdj <- -4
  layout(rbind(1:3,
               c(4:5, 0),
               rep(0, 3),
               6:8,
               9:11,
               rep(0, 3),
               12:14,
               15:17,
               18:20),
         heights = c(1, 1, 1/5, 1, 1, 1/5, 1, 1, 1))
  mlCF <- rowSums(concord) + postProb + iqStat[, "ufb"]
  Histy(concord[, "cluster"], cf = mlCF)
  Panel("a)", 0, yAdj)
  Histy(concord[, "mutual"], cf = mlCF)
  Histy(concord[, "quartet"], cf = mlCF)
  Histy(postProb, cf = mlCF, even = "log", breaks = 24)
  Histy(iqStat[, "ufb"], cf = mlCF, even = "log")
  
  iqCF <- rowSums(concord) + rowSums(iqStat)
  Histy(concord[, "cluster"], cf = iqCF)
  Panel("b)", 0, yAdj)
  Histy(concord[, "mutual"], cf = iqCF)
  Histy(concord[, "quartet"], cf = iqCF)
  Histy(iqStat[, "lbp"], cf = iqCF)
  Histy(iqStat[, "abayes"], cf = iqCF)
  Histy(iqStat[, "alrt"], cf = iqCF)
  
  tntCF <- rowSums(concord) + rowSums(tntStat) + bremer
  Histy(concord[, "cluster"], cf = tntCF)
  Panel("c)", 0, yAdj)
  Histy(concord[, "mutual"], cf = tntCF)
  Histy(concord[, "quartet"], cf = tntCF)
  Histy(tntStat[, "jak"], cf = tntCF)
  Histy(tntStat[, "boot"], cf = tntCF)
  Histy(tntStat[, "symFq"], cf = tntCF)
  Histy(tntStat[, "symGC"], cf = tntCF)
  Histy(tntStat[, "pois"], cf = tntCF)
  Histy(bremer, cf = tntCF)
  
  dev.off()
}

# ============================================================
# Appendix: NID vs. support metrics (referee request)
#
# Plots the normalized clustering information distance (NID)
# between each inferred edge and the reference tree (x-axis)
# against each support metric (y-axis), with Spearman rho
# (computed on all observations) annotated per panel.
# Metrics are ordered as in Fig. 2. Points are coloured by
# whether the split is in the reference tree (blue) or not
# (red); trend line shows the binned median across 40
# equal-width NID bins.
# ============================================================

# Metrics in the same order as Fig. 2 (unique, first appearance)
.fig_a1_metrics <- list(
  list(values = allDat$cluster,      name = "Clustering concordance"),
  list(values = allDat$mutual,       name = "Mutual clustering concordance"),
  list(values = allDat$quartet,      name = "Quartet concordance"),
  list(values = postProb,            name = "Posterior probability"),
  list(values = iqStat[, "ufb"],     name = "Ultra-fast bootstrap"),
  list(values = iqStat[, "lbp"],     name = "Local bootstrap prob."),
  list(values = iqStat[, "abayes"],  name = "Approx. Bayes"),
  list(values = iqStat[, "alrt"],    name = "Approx. LRT"),
  list(values = tntStat[, "jak"],    name = "TNT jackknife"),
  list(values = tntStat[, "boot"],   name = "TNT bootstrap"),
  list(values = tntStat[, "symFq"],  name = "Symmetric frequency"),
  list(values = tntStat[, "symGC"],  name = "Groups pres / cont"),
  list(values = tntStat[, "pois"],   name = "Poisson resampling"),
  list(values = bremer,              name = "Bremer support")
)

.nid <- 1 - partQual  # 0 = true split; > 0 = incorrect split

# Plot one panel: NID (x) vs. a support metric (y)
.NidPanel <- function(values, name, col_true = 3, col_false = 2,
                      n_sample = 2000) {
  ok      <- !is.na(values) & !is.na(.nid)
  x       <- .nid[ok]
  y       <- values[ok]
  correct <- partCorrect[ok]

  # Spearman correlation on all observations
  rho <- cor(x, y, method = "spearman")

  # Stratified sample for scatter display
  idx_t   <- which(correct)
  idx_f   <- which(!correct)
  keep    <- c(sample(idx_t, min(n_sample, length(idx_t))),
               sample(idx_f, min(n_sample, length(idx_f))))
  col_pts <- adjustcolor(ifelse(correct[keep], col_true, col_false),
                         alpha.f = 0.25)

  plot(x[keep], y[keep],
       pch        = 16,
       cex        = 0.3,
       col        = col_pts,
       xlim       = range(x, na.rm = TRUE),
       xlab       = "",
       ylab       = "",
       main       = name,
       frame.plot = FALSE)

  fit_gam <- gam(y ~ s(x, bs = "cs"),
                 family = gaussian(link = "log"),
                 method = "REML"
                 )
  plot_x  <- seq(0, max(x), length.out = 200)
  plot_y  <- predict(fit_gam, newdata = data.frame(x = plot_x),
                     type = "response")
  lines(plot_x, plot_y, lwd = 2)
  
  # Spearman rho annotation (top-right, inside panel)
  usr <- par("usr")
  text(usr[2], usr[4],
       labels = bquote(rho == .(sprintf("%.2f", rho))),
       adj    = c(1.05, 1.4),
       cex    = 0.7)
}

set.seed(4917)
cairo_pdf("Fig 3 - CID vs support.pdf", width = 7, height = 9)

layout(matrix(1:16, nrow = 4, ncol = 4, byrow = TRUE))
par(mar      = c(2.5, 2.5, 2, 0.5),
    oma      = c(2,   2,   0, 0),
    font.main = 1,
    cex.main  = 0.85,
    cex.axis  = 0.7,
    tcl       = -0.3,
    mgp       = c(2, 0.4, 0))

for (.m in .fig_a1_metrics) .NidPanel(.m$values, .m$name)

# Shared axis labels
mtext("Normalized clustering information distance",
      side = 1, outer = TRUE, line = 0.5, cex = 0.8)
mtext("Support value",
      side = 2, outer = TRUE, line = 0.5, cex = 0.8)

# Legend in the empty 15th slot
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend  = c("In reference tree", "Not in reference tree"),
       pch     = 16,
       col     = 3:2,
       bty     = "n",
       cex     = 0.9,
       pt.cex  = 1.2)

dev.off()
