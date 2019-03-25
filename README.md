# tafka-crisprQTL
Code associated with tafka-crisprQTL (Gasperini, et al. 2019)


## Scripts
The following scripts/files were used in differential expression testing. We also include information where to access relevant data files in the `Data` section below.

- `custom_monocle_da_functions.R`: other scripts source this very slightly modified version of monocle's differential expression testing functions (they report the intercept and the coefficient for the target term in the model). As such, this file must remain in the same directory as the scripts.

The scripts used for differential expression testing in the pilot and at-scale experiments are very similar but we keep them distinct for clarity:
- `get_deg.pilot.R`: differential expression script used for the pilot experiment

- `get_deg.at_scale.R`: differential expression script used for the at-scale experiment 

Each of these two scripts has the same interface, for example::
```
Rscript get_deg.pilot.R cds_object.rds \
grna_group \
gene_gRNAgroup_pair_table.pilot.txt \
output_file.txt
```

Note that for the `at-scale` experiment, in order to cut down on computational costs, we used a random subset of 50,000 cells from the experiment as the reference set rather than including all cells. This means that all cells with a guide to the specified target are included, but only cells in the 50,000 cell reference are included if they do not have a guide. This means that the call to this script is slightly different:
```
Rscript get_deg.pilot.R cds_object.rds \
grna_group \
gene_gRNAgroup_pair_table.at_scale.txt \
output_file.txt \
--reference_cells_rds 50k_reference_cells.rds
```

These scripts run the relevant DEG tests for a single target at a time (corresponding to `grna_group`). We then parallelized these over all grna_groups via our cluster. A list of all ggrna_group values for the pilot and at-scale screens are provided in the `Data` section below as a reference.

## Data
- `grna_groups.pilot.txt`: TODO desc
- `grna_groups.at_scale.txt`: TODO
- `50k_reference_cells.rds`: TODO
- `gene_gRNAgroup_pair_table.pilot.txt`: TODO
- `gene_gRNAgroup_pair_table.at_scale.txt`: TODO
- `all_deg_results.pilot.txt`: TODO
- `all_deg_results.at_scale`: TODO

## Required R Packages
The following packages are required to run the scripts:
```
library(argparse)
library(dplyr)
library(Matrix)
library(monocle)
library(stringr)
library(tidyr)
```

Note that we used Monocle 2. Because 1) monocle differential expression testing has not changed for some time and 2) we load our own differential expression functions that are modified to report intercepts and the target coefficient, monocle versioning should not have practical implications as long as the objects stored in our provided RDS files load without error.
