#!/usr/bin/env python
import os.path
import pandas as pd

configfile: "configuration/config.yaml"
sample_sheet = pd.read_csv(config["sample_sheet"],
    dtype=str,
    names = ["sample", "asm"]).set_index("sample")

wildcard_constraints:
    sample = "|".join(sample_sheet.index)

include: "rules/utils.smk"
include: "rules/1.annotate_BUSCOs.smk"
include: "rules/2.organize_sequences.smk"
include: "rules/3.reconstruct_gene_trees.smk"
include: "rules/4.reconstruct_species_tree.smk"

rule all:
    input:
        os.path.join(config["outdir"],"merged_best_single_copy_BUSCOs","summary.txt"),
        os.path.join(config["outdir"],"species_tree", config["project"] + ".speciesTree.txt")
    shell:
        """
        echo "Job done!"
        """