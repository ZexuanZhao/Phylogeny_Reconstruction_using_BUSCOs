rule separate_best_single_copy_BUSCOs:
    conda:
        os.path.join(workflow.basedir,"envs/biopython.yaml")
    input:
        full_table = os.path.join(config["outdir"],"compleasm","{sample}", config["busco_lineage"], "full_table.tsv"),
        proteins = os.path.join(config["outdir"],"compleasm","{sample}", config["busco_lineage"], "translated_protein.fasta")
    output:
        os.path.join(config["outdir"],"compleasm","{sample}", config["busco_lineage"], "best_single_copy_BUSCOs", "done.txt")
    params:
        script = "scripts/extract_and_rename_sequences.py",
        sample = "{sample}",
        outdir = os.path.join(config["outdir"],"compleasm","{sample}",config["busco_lineage"], "best_single_copy_BUSCOs")
    log:
        os.path.join(config["outdir"],"logs", "separate_best_single_copy_BUSCOs", "{sample}"+".log")
    threads:
        1
    shell:
        """
            python3 \
                {params.script} \
                {input.full_table} \
                {input.proteins} \
                {params.sample} \
                {params.outdir} > {log} 2>{log} && \
                touch {output} 
        """

rule merge_protein_files:
    conda:
        os.path.join(workflow.basedir,"envs/biopython.yaml")
    input:
        expand(os.path.join(config["outdir"],"compleasm","{sample}", config["busco_lineage"], "best_single_copy_BUSCOs", "done.txt"), sample = sample_sheet.index)
    output:
        os.path.join(config["outdir"], "merged_best_single_copy_BUSCOs", "summary.txt")
    params:
        script = "scripts/merge_faa_files.py",
        indirs = expand(os.path.join(config["outdir"],"compleasm","{sample}", config["busco_lineage"], "best_single_copy_BUSCOs"), sample = sample_sheet.index),
        outdir = os.path.join(config["outdir"], "merged_best_single_copy_BUSCOs"),
        n = config["minimum_sample_n"]
    log:
        os.path.join(config["outdir"],"logs", "merge_faa_files.log")
    threads:
        config["cpus"]
    shell:
        """
            python3 \
                {params.script} \
                {params.indirs} \
                {params.outdir} \
                -n {params.n} \
                -s {output} \
                -t {threads} \
                > {log} 2>{log} 
        """