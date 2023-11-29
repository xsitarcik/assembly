
rule blast__download_database:
    output:
        blast_db=protected(
            multiext(
                "{blast_db_dir}/{reference_tag}.{type}",
                "db",
                "hr",
                "in",
                "ot",
                "sq",
                "tf",
                "to",
            )
        ),
    params:
        blast_db_dir=lambda wildcards, output: os.path.dirname(output.blast_db[0]),
    log:
        "{blast_db_dir}/logs/{reference_tag}_{type}.log",
    retries: 3
    conda:
        "../envs/blast.yaml"
    shell:
        "(mkdir -p {params.blast_db_dir} && cd {params.blast_db_dir} && update_blastdb.pl --decompress {wildcards.reference_tag} && blastdbcheck -db {wildcards.reference_tag} ) > {log} 2>&1"


rule blast__create_header:
    output:
        header=temp("results/blast/header.tsv"),
    params:
        header="\t".join(BLAST_HEADER.split()),
    log:
        "logs/blast/header.log",
    conda:
        "../envs/echo.yaml"
    shell:
        "echo {params.header}  > {output.header}"


rule blast__query:
    input:
        blast_db=infer_blast_db,
        contigs="results/assembly/{sample}/contigs.fasta",
        header="results/blast/header.tsv",
    output:
        tsv="results/blast/{sample}/{reference_tag}.tsv",
        tsv_headerless=temp("results/blast/{sample}/{reference_tag}.tsv.tmp"),
    params:
        blast_db_dir=lambda wildcards, output: os.path.dirname(output.blast_db[0]),
        binary=infer_blast_binary,
        blast_format="6 {header}".format(header=BLAST_HEADER),
        max_number_of_hits=infer_max_number_of_hits,
    threads: min(config["threads"]["blast"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_blast,
    log:
        "logs/blast/{sample}/query_{reference_tag}.log",
    conda:
        "../envs/blast.yaml"
    shell:
        "(export BLASTDB={params.blast_db_dir} && {params.binary} -db {wildcards.reference_tag} -query {input.contigs} -out {output.tsv_headerless}"
        " -outfmt {params.blast_format:q} -num_threads {threads} -max_target_seqs {params.max_number_of_hits}"
        " && cat {input.header} {output.tsv_headerless} > {output.tsv}"
        " ) > {log} 2>&1"
