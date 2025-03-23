rule mafft:
    conda:
        os.path.join(workflow.basedir,"envs/mafft.yaml")
    input:
        os.path.join(config["outdir"],"merged_best_single_copy_BUSCOs","done.txt")
    output:
        os.path.join(config["outdir"],"aligned_merged_best_single_copy_BUSCOs","done.txt")
    params:
        indir=os.path.join(config["outdir"],"merged_best_single_copy_BUSCOs"),
        outdir=os.path.join(config["outdir"],"aligned_merged_best_single_copy_BUSCOs"),
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

rule iqtree:
    conda:
        os.path.join(workflow.basedir,"envs/iqtree.yaml")
    input:
        os.path.join(config["outdir"],"aligned_merged_best_single_copy_BUSCOs","done.txt")
    output:
        os.path.join(config["outdir"],"gene_trees","gene_trees.txt")
    params:
        indir = os.path.join(config["outdir"],"aligned_merged_best_single_copy_BUSCOs"),
        outdir = os.path.join(config["outdir"],"gene_trees"),
        n_parallele_iqtree = config["n_parallele_iqtree"]
    log:
        os.path.join(config["outdir"],"logs","iqtree.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.indir} -type f -name '*_aligned.faa' -print0 | \
                xargs -0 -P {params.n_parallele_iqtree} \
                    -I {{}} \
                    bash -c \
                        'iqtree \
                            -s "{{}}" \
                            -m MFP' \
                            > {log} 2>{log}
            cat {params.indir}/*.treefile > {output}
        """