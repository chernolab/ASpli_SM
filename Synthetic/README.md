# Supplementary material for ASpli2 paper

## Synthetic dataset
### Folders
- BAMs     : BAM files of the synthetic dataset (292 Mb)
- SimuData : Parameters and results of the simulated synthetic dataset
- ASpliRun : ASpli2 pre-computed analysis of the synthetic dataset (includes html event report)
- others   : Analysis of LeafCutter, MAJIQ and rMATS
- paperFigs: Folder where figures produced by the `run-comparison.R` script will be stored

### Scripts


- `run_ASpli.R`: Script to run the analysis from scratch. Results will be stored 
in the ASpliRun folder (new results will overwrite pre-calculated ones)
 
- `run-comparison.R`: Script to produce figures and tables of the reproducibility 
analysis presented in the paper from ASpli2 (*ASpliRun* folder) and 
other algorithms (*others* folder) pre-computed results.


