# Phylogeny Reconstruction using BUSCOs

## Description:
 - A snakemake pipeline to reconstruct phylogeny based on genome assemblies.
 - Annotate BUSCOs, sort and align sequences, and build gene trees and a species tree.

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 2 columns (no column name):
   - `sample` as the tip labels in the final species tree, `path to the genome assembly`
 - Modify configuration file - `configuration/config.yaml`:
   - `project`: a name for your project as the prefix of the final gene tree file.
   - `sample_sheet`: path to the sample sheet prepared above
   - `outdir`: path to the output directory
   - `busco_lineage`: lineage to annotate BUSCOs, check [here](https://busco-data.ezlab.org/v5/data/lineages/).
     - Note: use full lineage name, like xxx_odb12
   - `n_parallele_mafft`: number of jobs to run mafft in parallel
   - `n_parallele_iqtree`: number of jobs to run iqtree in parallel
   - `cpus`: number of CPUs to use
   - `astral`: name of the astral excutable in the `\bin` folder

## Environment:
 - Make sure snakemake is installed in current environment.
   
## Usage:
`snakemake --use-conda --cores [ncpu]`
