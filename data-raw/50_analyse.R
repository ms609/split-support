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

for (i in cli::cli_progress_along(seq_len(5), "Analysing")) {
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
  iqTree <- read.tree(IQFile(sim, aln, ".treefile"))
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
  partBrem <- tntTree[[1]]$node.label[as.numeric(names(as.Splits(tntTree[[1]]))) - NTip(tntTree[[1]])]
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
    
    ## TEMPORARY
    # if (file.info(concCache)$mtime < "2025-12-09 15:00:00 GMT") {
    #   conc[, "cluster"] <- ClusteringConcordance(partitions, dataset,
    #                                              normalize = FALSE)
    #   conc <- cbind(conc[, 1:5], "clusterNorm" = ClusteringConcordance(
    #     partitions, dataset, normalize = TRUE))
    #   write.table(conc, concCache)
    # }
    ## END TEMPORARY
    
  } else {
    cAll <- ClusteringConcordance(partitions, dataset, normalize = FALSE,
                                  return = "all")
    bestSums <- rowSums(cAll["hBest", , ])
    .Rezero <- function(value, zero) {
      (value - zero) / (1 - zero)
    }
    
    conc <- cbind(
      quartet = QuartetConcordance(partitions, dataset),
      cluster = rowSums(cAll["mi", , ]) / bestSums, # = ClustConc(norm = FALSE)
      phylo = PhylogeneticConcordance(partitions, dataset),
      mutual = MutualClusteringConcordance(partitions, dataset),
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


partCorrect <- do.call(c, partCorrectList)
partQual <- do.call(c, partQualList)
postProb <- do.call(c, postProbList)
concord <- do.call(rbind, concordList)
splitH <- do.call(rbind, splitHList)
bremer <- do.call(c, bremerList)
tntStat <- do.call(rbind, tntStatList)
iqStat <- do.call(rbind, iqStatList)

colnames(tntStat) <- c("symFq", "symGC", "boot", "jak", "pois")
colnames(iqStat) <- c("alrt", "lbp", "abayes", "ufb") # .iqtree output file gives order

common <- rowSums(is.na(concord)) == 0 &
  rowSums(is.na(tntStat)) == 0 &
  !is.na(bremer) &
  rowSums(is.na(iqStat)) == 0

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

CIndex <- function(score, target) {
  fit <- survcomp::concordance.index(
    x = score,
    surv.time  = -target,
    surv.event = rep(1, length(target)),
    method = "noether",
    na.rm  = TRUE)
  # fit$c.index is in [0,1]; 0.5 = random
  est <- fit$c.index
  se  <- fit$se
  ci  <- est + c(-1, 1) * 1.96 * se
  list(estimate = est, se = se, ci95 = ci, details = fit)
}
SomersD <- function(score, target) {
  ret <- CIndex(score, target)[c("estimate", "ci95")]
  lapply(ret, function(x) 2 * (x - 0.5))
}

SomersD(dat$postProb, dat$partQual)$estimate
SomersD(dat$postProb, dat$partQual)$ci95
SomersD(dat$quartet, dat$partQual)$ci95
SomersD(dat$cluster, dat$partQual)$ci95
SomersD(dat$phylo, dat$partQual)$ci95



# Platt scaling: logistic regression of occurs ~ score.

# devtools::install_github("bhklab/survcomp") 

# How well does a measure predict whether a split is in the true tree?
# We set `cf` to include only splits for which data is available under
# both `var` and `cf`, to allow a straight comparison.
Histy <- function(var, breaks = 16, even = TRUE, cf = var) { # "Mosaic plot"
  entries <- !is.na(var) & !is.na(cf)
  outcomes <- factor(partCorrect[entries], levels = c("FALSE", "TRUE"),
                     ordered = TRUE)
  var <- var[entries]
  brks <- if (even) {
    quantile(var, seq(0, 1, length.out = breaks))
  } else {
    seq(min(var), max(var), length.out = breaks)
  }
  brks <- unique(brks)
  pattern <- if (max(brks) > 2) "%.0f" else "%.2f"
  binLabels <- sprintf(pattern, brks[-length(brks)])
  bins <- cut(var, breaks = brks)
  
  col <- c("TRUE" = "3", "FALSE" = "2")
  if (substr(as.character(match.call()[-1])[[1]], 1, 7) != "concord") {
    col <- adjustcolor(col, alpha.f = 0.5)
  }
  title <- as.character(match.call()[2])
  title <- switch(
    title,
    "postProb" = "Posterior probability",
    "concord[, \"cluster\"]" = "Clustering concordance",
    "concord[, \"quartet\"]" = "Quartet concordance",
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
    title)
  
  tab <- table(bins, outcomes)
  #tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  #dimnames(tab)[[2]] <- rep("", dim(tab)[[2]])
  spTab <- spineplot(tab, main = title, col = col,
                     axes = FALSE,
                     xaxlabels = "", yaxlabels = "", xlab = "", ylab = "",
                     border = NA)
  
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
  
  text(x = x,
       y = -0.06,
       labels = binLabels,
       srt = 90,
       adj = 1,
       xpd = NA,
       cex = 0.8)
  
  
  roc <- pROC::roc(predictor = var, response = as.numeric(outcomes),
                   quiet = TRUE)
  sD <- SomersD(var, partQual[entries])
  
  message("n = ", sum(entries), ": ", title)
  # mtext(bquote(
  #   ROC-AUC == .(sprintf("%.3f", roc$auc)) * ";" ~
  #   D == .(sprintf("%.3f", sD$estimate))
  # ), 3, line = -0.3, cex = 0.6)
  mtext(paste0("ROC-AUC = ", sprintf("%.2f", roc$auc), "; ",
               "D = ", sprintf("%.2f", sD$estimate)),
        3, cex = 0.6)
}

{
  cairo_pdf("../char-concord/Fig 3 - edge concordance.pdf", 5.4, 8.4)
  par(mar = c(1.6, 1, 3, 1), font.main = 1, cex.main = 0.9)
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
  # Histy(concord[, "clusterNorm"], cf = postProb) # rubbish
  Histy(concord[, "mutual"], cf = mlCF)
  Histy(concord[, "quartet"], cf = mlCF)
  Histy(postProb, cf = mlCF)
  Histy(iqStat[, "ufb"], cf = mlCF)
  # Histy(concord[, "shared"], cf = postProb)
  # Histy(concord[, "phylo"], cf = postProb)
  #Histy(splitH, cf = postProb)
  
  iqCF <- rowSums(concord) + rowSums(iqStat)
  Histy(concord[, "cluster"], cf = iqCF)
  Histy(concord[, "mutual"], cf = iqCF)
  Histy(concord[, "quartet"], cf = iqCF)
  Histy(iqStat[, "lbp"], cf = iqCF)
  Histy(iqStat[, "abayes"], cf = iqCF)
  Histy(iqStat[, "alrt"], cf = iqCF)
  
  tntCF <- rowSums(concord) + rowSums(tntStat) + bremer
  Histy(concord[, "cluster"], cf = tntCF)
  Histy(concord[, "mutual"], cf = tntCF)
  Histy(concord[, "quartet"], cf = tntCF)
  Histy(bremer, cf = tntCF)
  Histy(tntStat[, "boot"], cf = tntCF)
  Histy(tntStat[, "jak"], cf = tntCF)
  Histy(tntStat[, "pois"], cf = tntCF)
  Histy(tntStat[, "symFq"], cf = tntCF)
  Histy(tntStat[, "symGC"], cf = tntCF)
  
  dev.off()
}









allCF <- rowSums(concord) + postProb + rowSums(tntStat) + rowSums(iqStat) + bremer
par(mfrow = c(5, 3), mar = rep(2, 4))
Histy(postProb, cf = allCF)
Histy(bremer, cf = allCF)
Histy(tntStat[, "symFq"], cf = allCF)
Histy(tntStat[, "symGC"], cf = allCF)
Histy(tntStat[, "boot"], cf = allCF)
Histy(tntStat[, "jak"], cf = allCF)
Histy(tntStat[, "pois"], cf = allCF)

Histy(iqStat[, "ufb"], cf = allCF)
Histy(iqStat[, "lbp"], cf = allCF)
Histy(iqStat[, "alrt"], cf = allCF)
Histy(iqStat[, "abayes"], cf = allCF)

Histy(concord[, "cluster"], cf = allCF)
# Histy(concord[, "clusterNorm"], cf = allCF)
Histy(concord[, "quartet"], cf = allCF)
Histy(concord[, "mutual"], cf = allCF)
# Histy(concord[, "shared"], cf = allCF)
# Histy(concord[, "phylo"], cf = allCF)


# The lower the Brier score is for a set of predictions,
# the better the predictions are calibrated.
# mclust::BrierScore(cbind(1 - postProb, postProb), partCorrect)
