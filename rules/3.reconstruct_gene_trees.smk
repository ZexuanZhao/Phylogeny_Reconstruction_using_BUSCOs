rule mafft:
    conda:
        os.path.join(workflow.basedir,"envs/mafft.yaml")
    input:
        os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged","summary.txt")
    output:
        os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned","done.txt")
    params:
        indir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged"),
        outdir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned"),
        n_parallele_mafft=config["n_parallele_mafft"],
        cpu_per_job = int(int(config["cpus"])/int(config["n_parallele_mafft"]))
    log:
        os.path.join(config["outdir"],"logs","mafft.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.indir} -type f -name '*.faa' -print0 | \
                xargs -0 -P {params.n_parallele_mafft} \
                    -I {{}} \
                    bash -c \
                        'mafft \
                            --thread {params.cpu_per_job} \
                            --auto \
                            "{{}}" \
                            > {params.outdir}/"$(basename "{{}}" .faa)"_aligned.faa' \
                            2>{log}
            touch {output}
        """

rule gblocks:
    conda:
        os.path.join(workflow.basedir,"envs/gblocks.yaml")
    input:
        os.path.join(config["outdir"], "best_single_copy_BUSCOs_merged_aligned", "done.txt")
    output:
        os.path.join(config["outdir"], "best_single_copy_BUSCOs_merged_aligned_gblocks", "done.txt")
    params:
        script = "scripts/Gblocks.py",
        indir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned"),
        outdir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned_gblocks"),
    log:
        os.path.join(config["outdir"],"logs","gblocks.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.indir} -type f -name '*_aligned.faa' -print0 | \
                xargs -0 -P {threads} \
                    -I {{}} \
                    bash -c \
                        'python scripts/Gblocks.py "{{}}"' \
                            1> {log} 2>{log}
            mv {params.indir}/*_aligned.faa.fsa {params.outdir}
            mv {params.indir}/*_aligned.faa.fsa.htm {params.outdir}
            touch {output}
        """
rule filter_sequences_with_Ns:
    conda:
        os.path.join(workflow.basedir,"envs/R.yaml")
    input:
        os.path.join(config["outdir"], "best_single_copy_BUSCOs_merged_aligned_gblocks", "done.txt")
    output:
        os.path.join(config["outdir"], "best_single_copy_BUSCOs_merged_aligned_gblocks_filtered", "done.txt")
    params:
        script = "scripts/filter_sequences_with_Ns.R",
        indir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned_gblocks"),
        outdir=os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned_gblocks_filtered"),
        p=config["p_Ns"]
    log:
        os.path.join(config["outdir"],"logs","filter_sequences_with_Ns.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.indir} -type f -name '*_aligned.faa.fsa' -print0 | \
                xargs -0 -P {threads} \
                    -I {{}} \
                    bash -c \
                        'Rscript {params.script} {params.p} "{{}}" {params.outdir}' \
                            1> {log} 2>{log}
            touch {output}
        """
rule iqtree:
    conda:
        os.path.join(workflow.basedir,"envs/iqtree.yaml")
    input:
        os.path.join(config["outdir"], "best_single_copy_BUSCOs_merged_aligned_gblocks_filtered", "done.txt")
    output:
        os.path.join(config["outdir"],"gene_trees","gene_trees.txt")
    params:
        indir = os.path.join(config["outdir"],"best_single_copy_BUSCOs_merged_aligned_gblocks_filtered"),
        outdir = os.path.join(config["outdir"],"gene_trees"),
        n_parallele_iqtree = config["n_parallele_iqtree"],
        cpu_per_job= int(int(config["cpus"]) / int(config["n_parallele_iqtree"]))
    log:
        iqtree_out=os.path.join(config["outdir"],"logs","iqtree.out"),
        iqtree_err=os.path.join(config["outdir"],"logs","iqtree.err"),
    threads:
        config["cpus"]
    shell:
        """
            find {params.indir} -type f -name '*_aligned.faa.filterN.fasta' -print0 | \
                xargs -0 -P {params.n_parallele_iqtree} \
                    -I {{}} \
                    bash -c \
                        'iqtree \
                            -redo \
                            -T {params.cpu_per_job} \
                            -s "{{}}" \
                            -m MFP' \
                            > {log.iqtree_out} 2>{log.iqtree_err}
            [ -d "{params.outdir}" ] || mkdir -p "{params.outdir}" 
            cat {params.indir}/*.treefile > {output} 
            find {params.indir} -name '*_aligned.faa.filterN.fasta.*' | \
                xargs mv --target-directory={params.outdir}
        """