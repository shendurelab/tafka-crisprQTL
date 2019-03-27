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
Rscript get_deg.pilot.R pilot_highmoi_screen.cds.rds \
grna_group \
gene_gRNAgroup_pair_table.pilot.txt \
output_file.txt
```

Note that for the `at-scale` experiment, in order to cut down on computational costs, we used a random subset of 50,000 cells from the experiment as the reference set rather than including all cells. This means that all cells with a guide to the specified target are included, but only cells in the 50,000 cell reference are included if they do not have a guide. This means that the call to this script is slightly different:
```
Rscript get_deg.at_scale.R at_scale_screen.cds.rds \
grna_group \
gene_gRNAgroup_pair_table.at_scale.txt \
output_file.txt \
--reference_cells_rds 50k_reference_cells.rds
```

These scripts run the relevant DEG tests for a single target at a time (corresponding to `grna_group`). We then parallelized these over all grna_groups via our cluster. A list of all ggrna_group values for the pilot and at-scale screens are provided in the `Data` section below as a reference.

## Data
Available GEO GSE120861.
- `at_scale_screen.cds.rds`: Processed scRNAseq CellDataSet for the at-scale screen (with metadata such as gRNAs assigned to each cell) to input into `get_deg.at_scale.R`.
- `pilot_highmoi_screen.cds.rds`: Processed scRNAseq CellDataSet for the pilot high MOI screen to input into `get_deg.pilot.R`. 
- `grna_groups.pilot.txt`: gRNAgroup identifiers for the pilot screen and their corresponding spacer sequences. The first column is a comprehensive list of all values to be used in place of `grna_group` in the `get_deg.pilot.R` command above. Only one group should be passed to this script at a time.
- `grna_groups.at_scale.txt`: gRNAgroup identifiers for the at-scale screen at corresponding spacer sequences. The first column is a comprehensive list of all values to be used in place of `grna_group` in the `get_deg.at_scale.R` command above. Only one group should be passed to this script at a time.
- `50k_reference_cells.rds`: File storing the names of the 50,000 cells used used for differential expression testing for the at-scale screen (passed as `--reference_cells_rds` to `get_deg.at_scale.R`).
- `gene_gRNAgroup_pair_table.pilot.txt`: Table of pilot screen's gene-target (gRNAgroup) pairs considered for differential expression testing (with corresponding information columns about each target and candidate target gene).
- `gene_gRNAgroup_pair_table.at_scale.txt`: Table of at-scale screen's gene-target (gRNAgroup) pairs considered for differential expression testing (with corresponding information columns about each target and candidate target gene).
- `all_deg_results.pilot.txt`: All differential expression results from the pilot screen, along with 'empirical p-values\*' (only calculated for `gRNA_groups` associated with decreases in candidate target gene's expression).
- `all_deg_results.at_scale`: All differential expression results from the at-scale screen, along with 'empirical p-values\*' (only calculated for `gRNA_groups` associated with decreases in candidate target gene's expression).
- `zero_inflated_outlier_genes_to_exclude.atscale.txt`: We observed that genes with a very high mean expression value relative to what you would expect given the number of cells in which they are expressed tended to show up frequently as false positives in permutation tests. These may be outliers with respect to their dispersion given that monocle uses shrinkage to estimate smoothed dispersion values when parameterizing the regression model used for differential gene expression. To avoid this source of false positives, we ignored this small set of outliers in downstream analysis. This set was derived from the data stored in `at_scale_screen.cds.rds`.
- `zero_inflated_outlier_genes_to_exclude.pilot.txt`: Same as above, but a list of outlier genes we identified from the `pilot_highmoi_screen.cds.rds`.

\* An empirical P-value was defined for each gene-gRNAgroup pair test that decreased the candidate target gene's expression. The empirical P-value was defined as: [(the number of NTCs with a smaller P-value than that test’s raw P-value) + 1] divided by [the total number of NTCs tests + 1].

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
