# Clustering concordance - simulation evaluation

This repository contains the simulation studies employed in Smith (forthcoming).

Analyses and results are contained within the `data-raw` directory.

To reproduce the simulation workflow, execute the `.R` scripts in numerical sequence.

## R scripts

- `_config.R`: Sets up the analysis and defines utility functions.
  Edit this file to specify:
  - The size of the tree and the datasets to be simulated.
  - The location of your own executables of
    [TNT](https://www.lillo.org.ar/phylogeny/tnt/),
    [MrBayes](https://nbisweden.github.io/MrBayes/download.html) and
    [IQ-TREE](https://iqtree.github.io/#download).
  - Your HPC login credentials.

- `10_simluate.R`: Simulates alignments using six rate categories drawn from a discretized gamma distribution.
  Outputs:
    - `reference-gam.tre`: Reference topology used for simulation, in newick format;
    - `alignments/gam###.nex`: Simulated alignment in NEXUS format; `###` denotes replicate ID.

- `20_

## References

Smith, Martin R. (forthcoming). "Which characters support which clades? Exploring the distribution of phylogenetic signal using concordant information."

