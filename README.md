# Benchmark of Differential TF activity methods

## Repository structure

Wrappers and higher-level functions are located in the `Scripts` folder.
Figures and downstream analyses are in the `Figures` folder.

Each other folder represents a benchmark processed separately:

* **fullFrags**: benchmark on real datasets, using all fragments
* **nucFree**: benchmark on real datasets, using only nucleosome-free fragments
* **fullOnNucFreePeaks**: benchmark on real datasets, using counts of all fragments on nucleosome-free peaks
* **simulations**: benchmark on semi-simulated datasets

Each of these folders contains one folder per dataset, with a `peaks` and a `seq_files` subfolders 
containing, respectively, the merged peaks and the (aligned reads). Due to size limitations the reads are not included in the repository, but are available on request.

To re-run the benchmark, simply do, from one of these folders:

```
source("../Scripts/runMethods.R")
datasets <- getDatasets()
runAll(datasets)
```

Results can be benchmarked using:

```
source("../Scripts/runMethods.R")
interactors <- readRDS("../Scripts/allInteractors.rds")
res <- compileBenchmark(datasets, interactors=interactors)
```
