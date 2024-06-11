# Benchmark of Differential TF activity methods

This repository contains the code used in Gerbaldo, Sonder et al., for the benchmark of differential TF activity methods based on ATAC-seq data.

## Repository structure

* **fullFrags**: benchmark on real datasets, using all fragments.
* **nucFree**: benchmark on real datasets, using only nucleosome-free fragments.
* **fullOnNucFreePeaks**: benchmark on real datasets, using counts of all fragments on nucleosome-free peaks.
* **withArchetypes**: benchmark on real datasets, using motif archetypes.
* **Scripts**: wrappers and higher-level functions.
* **Figures**: code behind the paper's figures.
* **simulations**: code for the generation of the semi-simulated datasets, and the benchmarks performed on those.
* **singleCell**: code for the cell-level analysis.
* **misc**: misc scripts preparing data objects for other analyses/plotting.
* **TRAFTAC**: code for the analysis of the TRAFTAC dataset.

## Main benchmark folders

Each of the 3 main benchmark folders contains one folder per dataset, with a `peaks` and a `seq_files` subfolders 
containing, respectively, the merged peaks and the (aligned) reads/fragments.
Due to size limitations the reads are not included in the repository, but are available on request.

To re-run the benchmark, simply do, from one of these folders:

```
source("../Scripts/runMethods.R")
datasets <- getDatasets()
runAll(datasets)
```

To run only the top methods, use:

```
runAll(datasets, methods=getMethods(onlyTop=TRUE))
```

Results can be benchmarked using:

```
source("../Scripts/compileBenchmark.R")
interactors <- readRDS("../Scripts/allInteractors.rds")
res <- compileBenchmark(datasets, interactors=interactors)
```

The precompiled results are also available in the respective folders.
