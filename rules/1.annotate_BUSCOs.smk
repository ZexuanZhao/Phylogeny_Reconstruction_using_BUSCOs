rule download_lineage:
    conda:
        os.path.join(workflow.basedir,"envs/compleasm.yaml")
    output:
        os.path.join(config["outdir"], "mb_downloads", config["busco_lineage"], "refseq_db.faa.gz")
    params:
        busco_lineage = config["busco_lineage"],
        path = os.path.join(config["outdir"], "mb_downloads")
    log:
        os.path.join(config["outdir"],"logs", "compleasm", "download.log")
    threads:
        1
    shell:
        """
            compleasm \
                download \
                -L {params.path} \
                {params.busco_lineage} \
                > {log} 2>{log}
        """

rule compleasm:
    conda:
        os.path.join(workflow.basedir, "envs/compleasm.yaml")
    input:
        unpack(get_asm),
        busco_lineage_downloaded_file = os.path.join(config["outdir"], "mb_downloads", config["busco_lineage"], "refseq_db.faa.gz")
    output:
        full_table=os.path.join(config["outdir"],"compleasm","{sample}",config["busco_lineage"],"full_table.tsv"),
        proteins=os.path.join(config["outdir"],"compleasm","{sample}",config["busco_lineage"],"translated_protein.fasta")
    params:
        outdir = os.path.join(config["outdir"], "compleasm", "{sample}"),
        busco_lineage = config["busco_lineage"],
        busco_lineage_dir = os.path.join(config["outdir"],"mb_downloads")
    log:
        os.path.join(config["outdir"],"logs", "compleasm", "{sample}.log")
    threads:
        config["cpus_per_compleasm"]
    shell:
        """
        compleasm \
            run \
            -a {input.asm} \
            -o {params.outdir} \
            -t {threads} \
            -l {params.busco_lineage} \
            -L {params.busco_lineage_dir} \
            > {log} 2>{log}
        """