# split-support: Agent Memory

## Project Overview

R-based simulation study evaluating character concordance measures and split
support metrics in phylogenetic inference. Accompanies the paper "Which
characters support which clades? Exploring the distribution of phylogenetic
signal using concordant information" (Smith, forthcoming).

## Directory Structure

```
split-support/
├── data-raw/
│   ├── _config.R                 # Central config: parameters, paths, helper fns
│   ├── 10_simulate.R             # Simulate alignments from reference tree
│   ├── 20_MrBayes*.R             # Bayesian inference (local + HPC variants)
│   ├── 30_iqtree.R               # Maximum-likelihood inference (IQ-TREE)
│   ├── 40_tnt.R                  # Parsimony inference (TNT)
│   ├── 80_byEdge.R               # Edge-wise analysis → Fig 2, Fig A1
│   ├── 90_byChar.R               # Character-wise analysis → Fig 3
│   ├── reference-gam.tre         # True reference tree (48 tips)
│   ├── alignments/               # 1 000 simulated NEXUS alignments
│   ├── MrBayes/                  # Bayesian results
│   ├── iqtree/                   # IQ-TREE results
│   ├── tnt/                      # TNT results
│   ├── concordance/              # Cached ClusteringConcordance results
│   └── entropy/                  # Cached ClusteringEntropy results
└── AGENTS.md
```

## Simulation Parameters (_config.R)

| Parameter | Value |
|-----------|-------|
| `nAln`    | 1 000 replicates |
| `nTip`    | 48 taxa |
| `nChar`   | 96 characters |
| `nCats`   | 6 gamma rate categories |
| `sim`     | `"gam"` (identifier used in all file-naming helpers) |

### File-naming helpers (from `_config.R`)

- `DataFile(sim, aln)` → NEXUS alignment path
- `MBFile(sim, aln, ext)` → MrBayes output path
- `IQFile(sim, aln, ext)` → IQ-TREE output path
- `TNTFile(sim, aln, type)` → TNT output path
- `ConcFile(sim, aln)` → concordance cache path
- `EntropyFile(sim, aln)` → entropy cache path
- `PartQFile(sim, aln)` → partition quality cache path

## Session Variables (populated by 80_byEdge.R)

| Variable | Type | Description |
|----------|------|-------------|
| `partCorrect` | logical vector | Whether each partition is in the reference tree |
| `partQual` | numeric vector | Normalised MCI of each partition vs. reference splits; 1 = true split |
| `postProb` | numeric vector | MrBayes posterior probabilities |
| `concord` | matrix | Concordance metrics: `quartet`, `cluster`, `phylo`, `mutual`, `shared`, `clusterNorm` |
| `tntStat` | matrix | TNT support: `symFq`, `symGC`, `boot`, `jak`, `pois` |
| `iqStat` | matrix | IQ-TREE support: `alrt`, `lbp`, `abayes`, `ufb` |
| `bremer` | numeric vector | Bremer (decay) support from TNT |
| `splitH` | matrix | Clustering entropy per partition |
| `allDat` | data.frame | All of the above combined (354 638 rows × 19 cols); includes NAs for method-specific metrics |
| `dat` | data.frame | `allDat` filtered to rows with no NAs in concordance + postProb columns |
| `refSplits` | Splits | Splits of the reference tree |
| `referenceTree` | phylo | The true 48-tip reference tree |

### Key derived quantity

```r
nid <- 1 - partQual  # Normalised clustering information distance (NID)
                     # 0 for every true split; (0, ~0.88] for incorrect splits
```

## Key Analyses (80_byEdge.R)

### Fig 2 — mosaic (spineplot) panels
Function `Histy(var, breaks, even, cf)` produces one mosaic panel showing how
`var` predicts `partCorrect`. Metrics are grouped by inference method (a =
MrBayes, b = IQ-TREE, c = TNT); concordance metrics appear in every group.
ROC-AUC and C-index are annotated on each panel.

**Fig 2 metric order** (unique, first appearance — also used for Fig A1):
`cluster`, `mutual`, `quartet`, `postProb`, `ufb`, `lbp`, `abayes`, `alrt`,
`jak`, `boot`, `symFq`, `symGC`, `pois`, `bremer`

### Fig A1 — NID vs. support scatter (appendix, referee request)
Function `.NidPanel(values, name)` plots NID (x) against a support metric (y).
Points are coloured by `partCorrect`; trend line is the binned median across 40
equal-width NID bins; Spearman ρ is annotated (computed on all observations,
not the display sample). Layout: 5 rows × 3 cols, 14 metric panels + 1 legend
slot. Output: `Fig A1 - NID vs support.pdf` (7 × 9 in).

## Concordance Metrics (TreeSearch package)

| Column | Function | Notes |
|--------|----------|-------|
| `cluster` | `ClusteringConcordance(norm=FALSE)` | Primary metric |
| `mutual` | `MutualClusteringConcordance()` | Primary metric |
| `quartet` | `QuartetConcordance()` | |
| `clusterNorm` | `ClusteringConcordance(norm=TRUE)` | Not used in study |
| `phylo` | `PhylogeneticConcordance()` | Not used in study |
| `shared` | `SharedPhylogeneticConcordance()` | Not used in study |

## Key Dependencies

**R packages:** `TreeTools`, `TreeSearch`, `TreeDist`, `phangorn`, `ape`,
`Hmisc` (C-index via `rcorr.cens`), `pROC`

**External software:** TNT, MrBayes 3.2.7, IQ-TREE 3.0.1

## Notes

- Concordance results are cached in `concordance/` and `entropy/` to avoid
  recomputation; delete cache files if partitions change.
- `partQual = 1` for **all** splits in the reference tree by construction
  (they are assigned 1 before the MCI calculation loop).
- The `common` logical vector flags rows where all methods have non-NA values
  (~37–45k rows); use this for fair cross-method comparisons.
- Systematic Biology figure guidelines require `frame.plot = FALSE` in base R
  plots (no surrounding box).
