# Load required libraries
library("TreeTools")
#library("TreeSearch")
devtools::load_all("../TreeSearch")

# Load configuration settings
source("data-raw/_config.R")

# Set simulation identifier here
sim <- "gam"

referenceTree <- file.path("data-raw", sprintf("reference-%s.tre", sim)) |>
  read.tree()
refSplits <- as.Splits(referenceTree)
tips <- names(read.nexus.data(DataFile(sim, "0001")))

# Eugh, I don't like growing vectors like this!
partCorrect <- logical(0)
postProb <- numeric(0)
concord <- numeric(0)
bremer <- numeric(0)
tntStats <- c("symFq", "symGC", "boot", "jak", "pois")
tntStat <- matrix(0, 0, length(tntStats), dimnames = list(NULL, tntStats))
ufb <- numeric(0)
iqStats <- c("alrt", "lbp", "abayes", "ufb") # .iqtree output file gives order
iqStat <- matrix(0, 0, length(iqStats), dimnames = list(NULL, iqStats))
splitH <- numeric(0)

for (i in cli::cli_progress_along(seq_len(nAln), "Analysing")) {
  aln <- alnIDs[[i]]
  
  # Load MrBayes partitions
  parts <- read.table(MBFile(sim, aln, "parts"), skip = 2 + nTip)
  partitions <- setNames(as.Splits(parts[, 2], tips), paste0("mb", parts[, 1]))
  
  # Load TNT partitions
  tntFile <- TNTFile(sim, aln, "ew")
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
    conc <- cbind(
      quartet = QuartetConcordance(partitions, dataset),
      cluster = ClusteringConcordance(partitions, dataset, normalize = FALSE),
      phylo = PhylogeneticConcordance(partitions, dataset),
      mutual = MutualClusteringConcordance(partitions, dataset),
      shared = SharedPhylogeneticConcordance(partitions, dataset),
      clusterNorm = ClusteringConcordance(partitions, dataset, normalize = TRUE)
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
  
  partCorrect <- c(partCorrect, partitions %in% refSplits)
  postProb <- c(postProb, pp, rep(0, sum(tntOnly, iqOnly, ufbOnly)))
  concord <- rbind(concord, conc)
  splitH <- rbind(splitH, h)
  bremer <- c(bremer, brem)
  tntStat <- rbind(tntStat, tntTags)
  iqStat <- rbind(iqStat, iqTags)
  
  stopifnot(dim(concord)[[1]] == dim(tntStat)[[1]])
  stopifnot(dim(concord)[[1]] == length(postProb))
}

common <- rowSums(is.na(concord)) == 0 &
  rowSums(is.na(tntStat)) == 0 &
  !is.na(bremer) &
  rowSums(is.na(iqStat)) == 0

model <- glm(partCorrect ~ postProb + concord + bremer + tntStat + iqStat,
             family = "binomial")

# How well does a measure predict whether a split is in the true tree?
# We set `cf` to include only splits for which data is available under
# both `var` and `cf`, to allow a straight comparison.
# We also report the R² and AIC of a binomial regression with split T/F.
Histy <- function(var, breaks = 20, even = TRUE, cf = var) { # "Mosaic plot"
  entries <- !is.na(var) & !is.na(cf)
  outcomes <- partCorrect[entries]
  var <- var[entries]
  if (even) {
    breaks <- quantile(var, seq(0, 1, length.out = breaks))
  }
  bins <- cut(var, breaks = unique(breaks))
  plot(
    table(bins, outcomes),
    main = as.character(match.call()[-1]),
    col = c("FALSE" = 2, "TRUE" = 3),
    xlab = "",
    ylab = ""
  )
  axis(1, signif(breaks), at = seq_along(breaks) / breaks, las = 2)
  
  # Predict whether a split is 'TRUE' using binomial regression
  m <- glm(outcomes ~ var, family = "binomial")
  smry <- summary(m)
  legend("left",
         text.font = 1,
         legend = c(paste("AIC:", round(smry$aic)),
                    paste("r2", signif(1 - (smry$deviance / smry$null.deviance), 3))))
}

par(mfrow = c(4, 2), mar = rep(2, 4))
Histy(postProb, cf = concord[, "quartet"])
Histy(concord[, "cluster"], cf = postProb)
Histy(concord[, "clusterNorm"], cf = postProb)
Histy(concord[, "quartet"], cf = postProb)
Histy(concord[, "mutual"], cf = postProb)
Histy(concord[, "shared"], cf = postProb)
Histy(concord[, "phylo"], cf = postProb)
Histy(splitH, cf = postProb)

par(mfrow = c(4, 3), mar = rep(2, 4))
Histy(bremer)
Histy(concord[, "quartet"], cf = bremer)
Histy(concord[, "mutual"], cf = bremer)
Histy(concord[, "shared"], cf = bremer)
Histy(concord[, "phylo"], cf = bremer)
Histy(concord[, "cluster"], cf = bremer)
Histy(concord[, "clusterNorm"], cf = bremer)

Histy(tntStat[, "symFq"])
Histy(tntStat[, "symGC"])
Histy(tntStat[, "boot"])
Histy(tntStat[, "jak"])
Histy(tntStat[, "pois"])

Histy(iqStat[, "ufb"])
Histy(iqStat[, "lbp"])
Histy(iqStat[, "alrt"])
Histy(iqStat[, "abayes"])

Peek <- function(var) {
  m <- glm(partCorrect ~ var, family = "binomial")
  smry <- summary(m)
  print(smry)
  print(paste("R2:", signif(1 - smry$deviance / smry$null.deviance)))
  AIC(m)
}
Peek(postProb)
Peek(bremer)
Peek(concord[, "quartet"])
Peek(concord[, "cluster"])
Peek(concord[, "clusterNorm"])
Peek(concord[, "mutual"])
Peek(iqStat[, "ufb"])

m <- glm(partCorrect ~ concord[, "quartet"] + postProb, family = "binomial")
AIC(m)
m <- glm(partCorrect ~ concord[, "quartet"], family = "binomial")

model <- glm(family = "binomial",
             partCorrect[common] ~ 
 #              postProb[common] +
#               bremer[common] +
               concord[common, "quartet"] +
               concord[common, "cluster"] +
               concord[common, "clusterNorm"] +
               concord[common, "phylo"] +
               concord[common, "mutual"] +
               concord[common, "shared"] +
               tntStat[common, "symFq"] +
               tntStat[common, "symGC"] +
               tntStat[common, "boot"] +
               tntStat[common, "jak"] + 
               tntStat[common, "pois"] +
               iqStat[common, "alrt"] +
               iqStat[common, "lbp"] +
               iqStat[common, "abayes"] +
               iqStat[common, "ufb"]
             )
step(model) # AIC
BIC(glm(family = "binomial", partCorrect[common] ~ postProb[common]))
BIC(glm(family = "binomial", partCorrect[common] ~ iqStat[common, "ufb"]))
BIC(glm(family = "binomial", partCorrect[common] ~ bremer[common]))
BIC(glm(family = "binomial", partCorrect[common] ~ tntStat[common, "symGC"]))
BIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "quartet"]))
BIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "mutual"]))
BIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "shared"]))
# BIC: https://stackoverflow.com/questions/19400494
step(model, criterion = "BIC", k = log(length(partCorrect)))

AIC(glm(family = "binomial", partCorrect[common] ~ postProb[common]))
AIC(glm(family = "binomial", partCorrect[common] ~ iqStat[common, "ufb"]))
AIC(glm(family = "binomial", partCorrect[common] ~ bremer[common]))
AIC(glm(family = "binomial", partCorrect[common] ~ tntStat[common, "symGC"]))
AIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "quartet"]))
AIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "mutual"]))
AIC(glm(family = "binomial", partCorrect[common] ~ concord[common, "shared"]))




# The lower the Brier score is for a set of predictions,
# the better the predictions are calibrated.
# mclust::BrierScore(cbind(1 - postProb, postProb), partCorrect)
