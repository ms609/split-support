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

# The two quartet measures reported in this study.  TreeSearch's defaults are
# now unit = "nrqs", chanceCorrect = TRUE, so every call pins all three
# arguments: a bare QuartetConcordance() would silently redefine the cached
# columns, and the cache guards below check only column names, not their
# meaning.
NaiveQC <- function(partitions, dataset, weight = TRUE) {
  QuartetConcordance(partitions, dataset, weight = weight,
                     unit = "quartet", chanceCorrect = FALSE)
}
NrqsQC <- function(partitions, dataset, weight = TRUE) {
  QuartetConcordance(partitions, dataset, weight = weight,
                     unit = "nrqs", chanceCorrect = TRUE)
}

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
      conc <- cbind(quartet = NaiveQC(partitions, dataset, weight = FALSE),
                    wQuartet = conc[, 1],

                    conc[, -1]
                    )
      write.table(conc, concCache)
    } else if (any(is.na(conc[, "quartet"]))) {
        conc[, "quartet"] <- NaiveQC(partitions, dataset, weight = FALSE)
        write.table(conc, concCache)
    }
    # Separate check, not an `else if`: a cache may need both migrations.
    if (!"wNrqs" %in% colnames(conc)) {
      conc <- cbind(conc, wNrqs = NrqsQC(partitions, dataset))
      write.table(conc, concCache)
    }
  } else {
    # For efficiency, calculate the complete concordance statistics once and
    # derive the associated measures below.
    # Output matches ClusteringConcordance(chanceCorrect = TRUE / FALSE)
    cAll <- ClusteringConcordance(
      partitions,
      dataset,
      chanceCorrect = FALSE,
      return = "all"
    )
    bestSums <- rowSums(cAll["hBest", , ])
    .Rezero <- function(value, zero) {
      (value - zero) / (1 - zero)
    }

    conc <- cbind(
      quartet = NaiveQC(partitions, dataset, weight = FALSE),
      wQuartet = NaiveQC(partitions, dataset),
      cluster = rowSums(cAll["mi", , ]) / bestSums, # = ClustConc(chanceC = F)
      phylo = PhylogeneticConcordance(partitions, dataset),
      mutual = MutualClusteringConcordance(partitions, dataset),
      # Not used in this study:
      shared = SharedPhylogeneticConcordance(partitions, dataset),
      clusterNorm = .Rezero(
        rowSums(cAll["mi", , ]) / bestSums,
        rowSums(cAll["miRand", , ]) / bestSums
      ), # = ClustConc(partitions, dataset, chanceCorrect = TRUE)
      # Appended last so that migrated and freshly written caches agree on
      # column order; all downstream access is by name in any case.
      wNrqs = NrqsQC(partitions, dataset)
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

here <- data.frame(
  conc = rowSums(is.na(concord)) == 0,
  tnt = rowSums(is.na(tntStat)) == 0,
  brem = !is.na(bremer),
  iq = rowSums(is.na(iqStat)) == 0)

hereFor <- paste0(
  ifelse(here[, "conc"], "C", ""),
  ifelse(here[, "tnt"], "T", ""),
  ifelse(here[, "brem"], "B", ""),
  ifelse(here[, "iq"], "I", ""))
table(hereFor)

herePal <- setNames(palette.colors(8),
                    c("", unique(hereFor)))

# Splits for which all four groups of statistics are available: concordance,
# TNT, Bremer and IQ-TREE -- i.e. hereFor == "CTBI". `here` has exactly those
# four logical columns, so a row sum of 4 selects them.
common <- rowSums(here) == 4


# Arrange in data.frame to allow subsequent filtering and analysis
allDat <- data.frame(
  occurs = partCorrect,
  partQual,
  postProb,
  tntStat,
  iqStat,
  bremer,
  concord,
  hereFor
)

dat <- data.frame(occurs = partCorrect, partQual, postProb, concord, hereFor) |> na.omit()

# Compute Somers' D, from which the C-index may be derived
SomersD <- function(score, target) {
  # survival::concordance() computes the same statistic as the
  # Hmisc::rcorr.cens() call this replaced -- verified identical to the last
  # digit on continuous, integer-valued, heavily tied and negative-valued
  # predictors against a partQual-shaped target -- but with an O(n log n)
  # balanced tree in place of rcorr.cens()'s O(n^2) double loop.  At
  # n = 354 638 that is 2 s against 19 min, and this function is called 51
  # times by the appendix figure alone.
  fit <- survival::concordance(target ~ score)

  # C index = (Dxy + 1) / 2, so Dxy = 2C - 1 and se(Dxy) = 2 se(C).
  est <- 2 * fit$concordance - 1
  # `$var` is an infinitesimal-jackknife variance.  rcorr.cens()'s "S.D." is
  # already a standard error of Dxy, not a standard deviation, and the two
  # agree to within 0.5%; the version of this function that divided "S.D." by
  # sqrt(n) understated the interval by a factor of sqrt(n) -- 595x at full
  # sample size.
  se <- 2 * sqrt(fit$var)

  ci  <- est + c(-1, 1) * 1.96 * se
  list(estimate = est, ci95 = ci, se = se)
}
# Derive the C-index from Somers' D.  `se` takes only the scale change, not the
# location shift, so the elements cannot be mapped over uniformly.
CIndex <- function(score, target) {
  d <- SomersD(score, target)
  list(estimate = (d$estimate + 1) / 2,
       ci95 = (d$ci95 + 1) / 2,
       se = d$se / 2)
}


# --- support chips: start ----------------------------------------------------
# Ranking 24 panels by eye is hard when the C-index and AUC are plain text, so
# each is printed on a rounded chip whose fill encodes the value.
#
# cividis, not viridis: it is monotonic in lightness, so "brighter is better"
# reads off directly, and its blue-to-yellow hue path stays clear of the
# red/green axis already spent on `partCorrect` in the mosaics.  It is also the
# CVD-optimised member of the family, which matters when hue carries meaning.
# The fills are fully opaque and the ink switches between black and white by
# luminance; a semi-transparent ramp would compress the dark end towards the
# page and cost the contrast that makes the small digits legible.
#
# Each statistic is scaled independently, anchored at 0.5 -- its
# no-discrimination value -- and stretched to its own attainable ceiling.  A
# shared [0.5, 1] domain would be tidier, but the C-index never exceeds 0.86
# here, so a fifth of the ramp would go unused and the 0.67-0.73 cluster would
# read as one flat grey.  Because the two domains differ, equal fills do NOT
# mean equal values across the two rows, and the legend must show both bars.
cIndexRange <- c(0.5, 0.9)
aucRange <- c(0.5, 1)
# Absolute character expansion wanted for the annotation: 0.5 of a 12 pt base,
# i.e. 6 pt, as the `mtext()` call this replaced produced.
chipCex <- 0.5
chipRamp <- viridisLite::cividis(256)

# `layout()` sets par("cex") to 0.66 for a grid this size, and text() /
# strwidth() multiply their `cex` by it whereas mtext() did not.  Divide it out,
# so the annotation keeps the size it has always had.
ChipCex <- function() chipCex / par("cex")

ChipFill <- function(x, range) {
  frac <- (x - range[[1]]) / diff(range)
  chipRamp[[1 + round(255 * min(max(frac, 0), 1))]]
}

# Whichever of black or white contrasts more strongly with the chip.  Contrast
# against the two is equal at a WCAG relative luminance of 0.179.
ChipInk <- function(fill) {
  channel <- col2rgb(fill)[, 1] / 255
  linear <- ifelse(channel <= 0.04045, channel / 12.92,
                   ((channel + 0.055) / 1.055) ^ 2.4)
  if (sum(linear * c(0.2126, 0.7152, 0.0722)) > 0.179) "#000000" else "#FFFFFF"
}

# Rounded rectangle.  Specified in inches, so the corner radius stays isotropic
# whatever user coordinates the panel happens to carry.
RoundRect <- function(x0, y0, x1, y1, radius, col, border = NA, n = 8) {
  radius <- min(radius, (x1 - x0) / 2, (y1 - y0) / 2)
  Corner <- function(cx, cy, from, to) {
    angle <- seq(from, to, length.out = n)
    cbind(cx + radius * cos(angle), cy + radius * sin(angle))
  }
  xy <- rbind(
    Corner(x1 - radius, y0 + radius, -pi / 2, 0),
    Corner(x1 - radius, y1 - radius, 0, pi / 2),
    Corner(x0 + radius, y1 - radius, pi / 2, pi),
    Corner(x0 + radius, y0 + radius, pi, 3 * pi / 2)
  )
  polygon(grconvertX(xy[, 1], "inches", "user"),
          grconvertY(xy[, 2], "inches", "user"),
          col = col, border = border, xpd = NA)
}

InchText <- function(x, y, labels, adj, col = "black", cex = ChipCex()) {
  text(grconvertX(x, "inches", "user"), grconvertY(y, "inches", "user"),
       labels = labels, adj = adj, col = col, cex = cex, xpd = NA)
}

# One annotation row: `prefix`, then `value` on its chip, then `suffix`.  `x`
# and `y` are in inches from the foot of the device; `xAdj` anchors the row
# horizontally (0 = left edge at `x`, 0.5 = centred, 1 = right edge at `x`).
# `prefix` may be a plotmath expression as well as a string.
ChipRow <- function(x, y, prefix, value, range, suffix = "", xAdj = 0.5,
                    cex = ChipCex()) {
  fill <- ChipFill(value, range)
  label <- sprintf("%.2f", value)
  digitW <- strwidth("0", units = "inches", cex = cex)
  digitH <- strheight("0", units = "inches", cex = cex)
  padX <- 0.5 * digitW
  padY <- 0.4 * digitH
  wPrefix <- strwidth(prefix, units = "inches", cex = cex)
  wSuffix <- strwidth(suffix, units = "inches", cex = cex)
  wChip <- strwidth(label, units = "inches", cex = cex) + 2 * padX
  chipL <- x - xAdj * (wPrefix + wChip + wSuffix) + wPrefix
  chipR <- chipL + wChip
  RoundRect(chipL, y - digitH / 2 - padY, chipR, y + digitH / 2 + padY,
            radius = 0.4 * (digitH + 2 * padY), col = fill)
  InchText(chipL, y, prefix, adj = c(1, 0.5), cex = cex)
  InchText((chipL + chipR) / 2, y, label, adj = c(0.5, 0.5), ChipInk(fill),
           cex = cex)
  if (is.character(suffix) && nzchar(suffix)) {
    InchText(chipR, y, suffix, adj = c(0, 0.5), cex = cex)
  }
  # Height of the chip, so callers can stack rows without re-deriving it
  invisible(digitH + 2 * padY)
}

# Height of a chip in inches, for laying out a stack before drawing it
ChipHeight <- function(cex = ChipCex()) {
  1.8 * strheight("0", units = "inches", cex = cex)
}

# The panel title is drawn here rather than by `spineplot()`, whose default
# placement leaves too little room below it for the chips.  Both statistics are
# passed as (lower, estimate, upper), the order `pROC::ci.auc()` returns.
SupportAnnotation <- function(mainTitle, cIndexCI, aucCI) {
  # par("mai") / par("mar") converts margin lines to inches exactly
  lineIn <- par("mai")[[3]] / par("mar")[[3]]
  top <- grconvertY(par("usr")[[4]], "user", "inches")
  centre <- grconvertX(mean(par("usr")[1:2]), "user", "inches")
  graphics::title(main = mainTitle, line = 2.6)
  # AUC first, to match the order the two statistics are introduced in the text;
  # `ChipLegend()` stacks its bars the same way.
  # "C-idx", not "C-index": the interval takes the row to the full width of the
  # panel, and the two rows must not disagree on how many decimals they show.
  ChipRow(centre, top + 1.4 * lineIn, "AUC = ", aucCI[[2]], aucRange,
          Bounds(aucCI))
  ChipRow(centre, top + 0.45 * lineIn, "C-idx = ", cIndexCI[[2]], cIndexRange,
          Bounds(cIndexCI))
}

# The 95% interval of a (lower, estimate, upper) triple, as a chip suffix.
# Three decimals: every interval here is narrower than the second decimal
# place, so two would print a bound equal to the estimate.
Bounds <- function(ci) sprintf(" (%.3f\U2013%.3f)", ci[[1]], ci[[3]])

# One bar of a colour key, in the user coordinates of a 0-1 window.  With
# `labelSide = "top"` the label sits above the bar and the bounds beneath its
# ends; with "left" the label sits to the left and the bounds run to the right,
# which is more compact where several bars must share a slot.
ChipBar <- function(x0, x1, y, height, range, label, cex = ChipCex(),
                    labelSide = "top") {
  rasterImage(t(chipRamp), x0, y, x1, y + height, interpolate = TRUE)
  rect(x0, y, x1, y + height, border = "#00000060", lwd = 0.5, xpd = NA)
  if (labelSide == "top") {
    text(x0, y + height + 0.02, label, adj = c(0, 0), cex = cex, xpd = NA)
    text(x0, y - 0.02, sprintf("%.2f", range[[1]]), adj = c(0, 1), cex = cex,
         xpd = NA)
    text(x1, y - 0.02, sprintf("%.2f", range[[2]]), adj = c(1, 1), cex = cex,
         xpd = NA)
  } else {
    text(x0 - 0.03, y + height / 2, label, adj = c(1, 0.5), cex = cex,
         xpd = NA)
    text(x1 + 0.03, y + height / 2,
         sprintf("%.2f\U2013%.2f", range[[1]], range[[2]]), adj = c(0, 0.5),
         cex = cex, xpd = NA)
  }
}

# Colour key for Fig 4's chips, drawn into a spare slot of the layout.  One bar
# per statistic, because the two are scaled to different ceilings.
ChipLegend <- function() {
  op <- par(mar = c(1, 0.5, 1, 0.5))
  on.exit(par(op))
  plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  ChipBar(0, 1, 0.66, 0.13, aucRange, "AUC")
  ChipBar(0, 1, 0.28, 0.13, cIndexRange, "C-index")
  arrows(0.72, 0.05, 1, 0.05, length = 0.02, lwd = 0.5, xpd = NA)
  text(0.68, 0.05, "better", adj = c(1, 0.5), cex = ChipCex(), xpd = NA)
}


# Fig 3's three C-indices, chipped on the same cividis ramp but each with its
# own bounds: they occupy very different parts of the [0, 1] interval (observed
# spans 0.62-0.86, 0.56-0.93, 0.38-0.67), and one shared domain would flatten
# all three.  Each domain is anchored on 0.5, the value a measure that ranks
# edges no better than chance would take:
#
#   C-index   [0.5, 0.9]  over all edges carrying data for the measure.  The
#             same domain as Fig 4's C-index chips, deliberately, so that a
#             chip means the same thing in both figures.
#   C_CTBI    [0.5, 0.95] over `common` edges only; a higher ceiling because
#             posterior probability reaches 0.93 on this subset.
#   C_CTBI'   [0.3, 0.7]  over incorrect `common` edges only.  This one is
#             *symmetric* about 0.5, because it is the one statistic that goes
#             below it -- the site concordance factor scores 0.38, i.e. ranks
#             incorrect edges in the wrong order.  Symmetry buys an absolute
#             reading: mid-grey is uninformative, brighter is informative, and
#             darker than mid-grey is actively inverted.
#
# Bounds are hard-coded rather than derived at draw time so the key reads the
# same from one run to the next; `.NidChips()` warns if a value falls outside,
# which `ChipFill()` would otherwise clamp silently to the end of the ramp.
nidCtbiRange <- c(0.5, 0.95)
nidCtbi1Range <- c(0.3, 0.7)
nidChipCex <- 0.78 # as the plain-text annotation this replaced used

# Three chipped rows down the top-right corner of a Fig 3 panel
.NidChips <- function(cIndex, ctbi, ctbi1) {
  usr <- par("usr")
  right <- grconvertX(usr[[2]], "user", "inches")
  top <- grconvertY(usr[[4]], "user", "inches")
  pitch <- 1.15 * ChipHeight(nidChipCex)
  inset <- 0.5 * strwidth("0", units = "inches", cex = nidChipCex)
  rows <- list(
    list(value = cIndex, range = cIndexRange, prefix = "C-index = "),
    list(value = ctbi, range = nidCtbiRange,
         prefix = as.expression(bquote(C[CTBI] == ""))),
    list(value = ctbi1, range = nidCtbi1Range,
         prefix = as.expression(bquote(C[CTBI * "'"] == "")))
  )
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    if (row$value < row$range[[1]] || row$value > row$range[[2]]) {
      warning("Chip value ", signif(row$value, 3), " outside [",
              row$range[[1]], ", ", row$range[[2]], "]; clamped")
    }
    ChipRow(right - inset, top - (i - 0.7) * pitch, row$prefix, row$value,
            row$range, xAdj = 1, cex = nidChipCex)
  }
}

# Colour key for Fig 3's chips: three bars, sharing the legend slot with the
# "Stats available" point key above them.
NidChipLegend <- function() {
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  bars <- list(
    list(range = cIndexRange, label = "C-index"),
    list(range = nidCtbiRange, label = as.expression(bquote(C[CTBI]))),
    list(range = nidCtbi1Range, label = as.expression(bquote(C[CTBI * "'"])))
  )
  for (i in seq_along(bars)) {
    ChipBar(0.42, 0.72, 0.23 - (i - 1) * 0.09, 0.05, bars[[i]]$range,
            bars[[i]]$label, cex = nidChipCex, labelSide = "left")
  }
}
# --- support chips: end ------------------------------------------------------


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
    "concord[, \"quartet\"]" = "Quartet concordance (unwtd)",
    "concord[, \"wQuartet\"]" = "Quartet concordance",
    "concord[, \"wNrqs\"]" = "NRQS concordance",
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
    main = "", # drawn by SupportAnnotation(), which needs the room below it
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

  cacheFile <- file.path(rocDir,
                         gsub("[ ,\"\\[]|\\]", "",
                              paste(call[-1], collapse = "-")))
  message(cacheFile)
  haveRoc <- FALSE
  haveCIdx <- FALSE
  if (file.exists(cacheFile)) {
    load(cacheFile)
    haveRoc <- TRUE
    # Caches written before the standard-error fix hold a `ci95` that is
    # sqrt(n) too narrow, and carry no `se` element.  Recompute those: it costs
    # seconds now that SomersD() is O(n log n), and leaving wrong intervals in
    # a cache is a trap.
    haveCIdx <- !is.null(cIdx$se)
  }
  if (!haveCIdx) {
    cIdx <- CIndex(var, partQual[entries])
  }
  if (!haveRoc) {
    roc <- pROC::roc(predictor = var, response = as.numeric(outcomes),
                     quiet = TRUE)
  }
  if (!haveRoc || !haveCIdx) {
    save(roc, cIdx, file = cacheFile)
  }

  message("n = ", sum(entries), ": ", title)

  rocVal <- as.numeric(pROC::ci.auc(roc))

  SupportAnnotation(title,
                    c(cIdx$ci95[[1]], cIdx$estimate, cIdx$ci95[[2]]),
                    rocVal)
}

# Produce figure as PDF
{
  cairo_pdf("Fig 4 - edge concordance.pdf", 5.4, 8.4)
  par(mar = c(1.6, 0.8, 4, 0.8), font.main = 1, cex.main = 0.9)
  yAdj <- -4
  layout(rbind(1:4,
               c(0, 5:6, 0),
               rep(0, 4),
               7:10,
               11:14,
               rep(0, 4),
               15:18,
               19:22,
               # Slot 25 is the chip colour key: bottom right, out of the way
               c(0, 23:24, 25)),
         heights = c(1, 1, 1/5, 1, 1, 1/5, 1, 1, 1))
  mlCF <- rowSums(concord) + postProb + iqStat[, "ufb"]
  Histy(concord[, "cluster"], cf = mlCF)
  Panel("a)", 0, yAdj)
  Histy(concord[, "mutual"], cf = mlCF)
  Histy(concord[, "wNrqs"], cf = mlCF)
  Histy(concord[, "wQuartet"], cf = mlCF)
  Histy(postProb, cf = mlCF, even = "log", breaks = 24)
  Histy(iqStat[, "ufb"], cf = mlCF, even = "log")
  
  iqCF <- rowSums(concord) + rowSums(iqStat)
  Histy(concord[, "cluster"], cf = iqCF)
  Panel("b)", 0, yAdj)
  Histy(concord[, "mutual"], cf = iqCF)
  Histy(concord[, "wNrqs"], cf = iqCF)
  Histy(concord[, "wQuartet"], cf = iqCF)
  Histy(iqStat[, "lbp"], cf = iqCF)
  Histy(iqStat[, "abayes"], cf = iqCF)
  Histy(iqStat[, "alrt"], cf = iqCF)
  Histy(iqStat[, "sCF"], cf = iqCF)
  
  tntCF <- rowSums(concord) + rowSums(tntStat) + bremer
  Histy(concord[, "cluster"], cf = tntCF)
  Panel("c)", 0, yAdj)
  Histy(concord[, "mutual"], cf = tntCF)
  Histy(concord[, "wNrqs"], cf = tntCF)
  Histy(concord[, "wQuartet"], cf = tntCF)
  Histy(tntStat[, "jak"], cf = tntCF)
  Histy(tntStat[, "boot"], cf = tntCF)
  Histy(tntStat[, "symFq"], cf = tntCF)
  Histy(tntStat[, "symGC"], cf = tntCF)
  Histy(tntStat[, "pois"], cf = tntCF)
  Histy(bremer, cf = tntCF)
  ChipLegend()

  dev.off()
}

# ============================================================
# Fig 3: NID vs. support metrics (referee request)
#
# Plots the normalized clustering information distance (NID)
# between each inferred edge and the reference tree (x-axis)
# against each support metric (y-axis), with three C-indices
# chipped in the top-right corner of each panel.  Metrics are
# ordered as in Fig. 4.  Points are coloured by which methods
# carry data for the edge; the two trend lines are GAM fits,
# over `common` edges and over all edges with data.
# ============================================================

# Metrics in the same order as Fig. 4 (unique, first appearance)
.fig_a1_metrics <- list(
  list(values = allDat$cluster,      name = "Clustering concordance"),
  list(values = allDat$clusterNorm,  name = "Adjusted clustering conc."),
  list(values = allDat$mutual,       name = "Mutual clustering concordance"),
  list(values = allDat$wNrqs,        name = "NRQS concordance"),
  list(values = allDat$wQuartet,     name = "Quartet concordance"),
  list(values = postProb,            name = "Posterior probability"),
  list(values = iqStat[, "ufb"],     name = "Ultra-fast bootstrap"),
  list(values = iqStat[, "lbp"],     name = "Local bootstrap prob."),
  list(values = iqStat[, "abayes"],  name = "Approx. Bayes"),
  list(values = iqStat[, "alrt"],    name = "Approx. LRT"),
  list(values = iqStat[, "sCF"],     name = "Site concordance factor"),
  list(values = tntStat[, "jak"],    name = "TNT jackknife"),
  list(values = tntStat[, "boot"],   name = "TNT bootstrap"),
  list(values = tntStat[, "symFq"],  name = "Symmetric frequency"),
  list(values = tntStat[, "symGC"],  name = "Groups pres / cont"),
  list(values = tntStat[, "pois"],   name = "Poisson resampling"),
  list(values = bremer,              name = "Bremer support")
)

.nid <- 1 - partQual  # 0 = true split; > 0 = incorrect split

# Everything expensive in a panel -- three C-indices and two GAM fits -- is
# cached per metric, as `Histy()` caches its ROC and C-index above.  Only the
# scatter sample is redrawn, so a second run of this figure costs seconds
# rather than the ~2 h it took when the C-index was quadratic and nothing was
# stored.
#
# The guard checks the file's version and nothing else, so bump
# `nidCacheVersion` after ANY change to what is computed below: the GAM family,
# `gamma`, the spline basis, the subsets the C-indices use, or the meaning of a
# metric column.  Deleting `data-raw/nid/` has the same effect.
nidCacheVersion <- 1L

.NidCacheFile <- function(name) {
  file.path(nidDir, gsub("[^A-Za-z0-9]", "", name))
}

# Two metrics whose names differed only in punctuation would silently share a
# cache file, and the version guard would not notice.
stopifnot(!anyDuplicated(
  vapply(.fig_a1_metrics, function(m) .NidCacheFile(m$name), "")
))

.NidCache <- function(name, Compute) {
  cacheFile <- .NidCacheFile(name)
  if (file.exists(cacheFile)) {
    cached <- readRDS(cacheFile)
    if (identical(cached$version, nidCacheVersion)) {
      return(cached)
    }
  }
  message("Computing ", name)
  cached <- c(Compute(), list(version = nidCacheVersion))
  saveRDS(cached, cacheFile)
  cached
}

# Plot one panel: NID (x) vs. a support metric (y)
.NidPanel <- function(values, name, n_sample = 2000) {
  ok      <- !is.na(values) & !is.na(.nid)
  x       <- .nid[ok]
  y       <- values[ok]
  correct <- partCorrect[ok]
  col     <- herePal[hereFor[ok]]

  stats <- .NidCache(name, function() {
    .CalcCI <- function(subset) CIndex(values[subset], partQual[subset])
    # `x` and `y` are already subset by `ok`, so `common` (full length) must be
    # subset the same way before it can index them.
    okCommon <- common[ok]
    xc <- x[okCommon]
    plotX <- seq(0, max(x), length.out = 200)
    gamma <- 2

    allGam <- mgcv::gam(y[okCommon] ~ s(xc, bs = "cs"),
                        family = gaussian(link = "log"),
                        gamma = gamma,
                        method = "REML"
                        )
    fitGam <- mgcv::gam(y ~ s(x, bs = "cs"),
                        family = gaussian(link = "log"),
                        gamma = gamma,
                        method = "REML"
                        )
    list(
      cidx = .CalcCI(ok),
      all_ci = .CalcCI(common),
      all_ci_1 = .CalcCI(common & !partCorrect),
      plotX = plotX,
      allY = predict(allGam, newdata = data.frame(xc = plotX),
                     type = "response"),
      plotY = predict(fitGam, newdata = data.frame(x = plotX),
                      type = "response")
    )
  })

  # Stratified sample for scatter display
  idx_t   <- which(correct)
  idx_f   <- which(!correct)
  keep    <- c(sample(idx_t, min(n_sample, length(idx_t))),
               sample(idx_f, min(n_sample, length(idx_f))))
  col_pts <- adjustcolor(col, alpha.f = 0.25)

  plot(x[keep], y[keep],
       pch        = 16,
       cex        = 0.4,
       col        = col_pts,
       xlim       = range(x, na.rm = TRUE),
       xlab       = "",
       ylab       = "",
       main       = name,
       frame.plot = FALSE)

  lines(stats$plotX, stats$allY, lwd = 2, col = herePal["CTBI"])
  lines(stats$plotX, stats$plotY, lwd = 2)

  # C-index annotations, chipped, down the top-right corner inside the panel
  .NidChips(stats$cidx$estimate, stats$all_ci$estimate,
            stats$all_ci_1$estimate)
}

set.seed(4917)
cairo_pdf("Fig 3 - CID vs support.pdf", width = 7, height = 9)

layout(matrix(1:18, nrow = 6, ncol = 3, byrow = TRUE))
par(mar      = c(2.5, 2.5, 2, 0.5),
    oma      = c(2,   2,   0, 0),
    font.main = 1,
    cex.main  = 0.85,
    cex.axis  = 0.7,
    tcl       = -0.3,
    mgp       = c(2, 0.4, 0))

for (.m in .fig_a1_metrics) .NidPanel(.m$values, .m$name)

# Shared axis labels
mtext("Edge support value",
      side = 2, outer = TRUE, line = 0.5, cex = 0.8)
mtext("Normalized clustering information distance to true tree",
      side = 1, outer = TRUE, line = 0.5, cex = 0.8)

hereIds <- names(herePal[-1])[c(3, 2, 5, 4, 6, 7, 1)]
# Legend in the empty slot.  The point key is pushed to the top and shrunk a
# little to leave room for the chip colour bars beneath it.
par(mar = c(0, 0, 0, 0))
plot.new()
legend("top",
       title = "Stats available:",
       legend  = paste0(hereIds, " (", table(hereFor)[hereIds], ")"),
       pch     = 16,
       col     = herePal[hereIds],
       bty     = "n",
       cex     = 0.9,
       pt.cex  = 1.2)
NidChipLegend()

dev.off()
