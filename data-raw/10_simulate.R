source("data-raw/_config.R")

# Simulate all characters with identical rates
# Create reference tree
set.seed(0)
referenceTree <- ape::rtree(nTip, equiprob = TRUE)
referenceTree$tip.label <- paste0("t", seq_len(nTip))
referenceTree <- TreeTools::RootTree(referenceTree, "t1")
treeLength <- sum(referenceTree$edge.length)
rate <- 12 / treeLength
print(signif(rate)) # mb.nex: prset brlenspr=unconstrained:uniform(0,<RATE>);
if (interactive()) plot(referenceTree)
ape::write.tree(referenceTree, file = "data-raw/reference-aln.tre")

# Simulate J-C alignments
for (i in formatC(1:1000, width = 4, flag = "0")) {
  ape::write.nexus.data(
    toupper(TreeTools::PhyDatToMatrix(
      phangorn::simSeq(referenceTree, nChar,
                       rootseq = rep("a", nChar),
                       rate = rate
                       ) # Jukes-Cantor model
    )),
    file = paste0("data-raw/alignments/aln", i, ".nex"),
    interleaved = FALSE 
  )
}

# Simulate characters under six gamma-distributed rate categories
# Create reference tree
set.seed(1984)
referenceTree <- ape::rtree(nTip, equiprob = TRUE)
referenceTree$tip.label <- paste0("t", seq_len(nTip))
referenceTree <- TreeTools::RootTree(referenceTree, "t1")
treeLength <- sum(referenceTree$edge.length)
rate <- 12 / treeLength
print(signif(rate)) # mb.nex: prset brlenspr=unconstrained:uniform(0,<RATE>);
if (interactive()) plot(referenceTree)
ape::write.tree(referenceTree, file = "data-raw/reference-gam.tre")

if (nChar %% nCats != 0) {
  warning("Remainder alert: ",
          "can't divide characters evenly between rate categories")
}
charPerCat <- nChar / nCats
cats <- phangorn::discrete.gamma(1, nCats)

set.seed(1984) # Reset seed for character simulation
# Simulate J-C alignments
for (i in formatC(1:1000, width = 4, flag = "0")) {
  do.call(cbind, lapply(cats, function(cat) {
    phangorn::simSeq(# Jukes-Cantor model
      referenceTree,
      charPerCat,
      rootseq = rep("a", charPerCat),
      rate = cat * rate
    ) |>
      TreeTools::PhyDatToMatrix()
    })) |> toupper() |>
    ape::write.nexus.data(file = paste0("data-raw/alignments/gam", i, ".nex"),
                     interleaved = FALSE)
}



# TODO: Consider simulating data under a non-JC model so that MCMC analysis
# uses a mis-specified model.

# See: http://www.iqtree.org/doc/AliSim
