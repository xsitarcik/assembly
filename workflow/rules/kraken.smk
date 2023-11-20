rule curl__download_kraken_db:
    output:
        protected("{prefix_dir}/hash.k2d"),
    params:
        url=lambda wildcards, output: "https://genome-idx.s3.amazonaws.com/kraken/{tag}.tar.gz".format(
            tag=os.path.basename(os.path.dirname(output[0]))
        ),
        dirpath=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
    log:
        "{prefix_dir}/logs/download.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "(mkdir -p {params.dirpath} && curl -SL {params.url} | tar zxvf - -C {params.dirpath}) > {log} 2>&1"


rule kraken__analysis:
    input:
        kraken_tax=os.path.join(config["kraken_dir"], "hash.k2d"),
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
    output:
        kraken_output=temp("results/kraken/{sample}.kraken"),
        report="results/kraken/{sample}.kreport2",
    params:
        save_memory="--memory-mapping" if config["kraken__params"]["save_memory"] else "",
        db_dir=lambda wildcards, input: os.path.dirname(input.kraken_tax),
    threads: min(config["threads"]["kraken"], config["max_threads"])
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "(kraken2 --db {params.db_dir} --threads {threads} --paired --gzip-compressed"
        " {params.save_memory} --report {output.report} {input.r1} {input.r2} 1> {output.kraken_output}) 2> {log}"


rule krona__update_taxonomy:
    output:
        protected("{prefix_dir}/taxonomy.tab"),
    params:
        tax_dir=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
    log:
        "{prefix_dir}/logs/update_taxonomy.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "ktUpdateTaxonomy.sh {params.tax_dir} 1> {log} 2>&1"


rule kraken__krona_chart:
    input:
        kraken_output="results/kraken/{sample}.kreport2",
        tax_tab=os.path.join(config["krona_dir"], "taxonomy.tab"),
    output:
        report(
            "results/kraken/kronas/{sample}.html",
            category="{sample}",
            labels={
                "Type": "Krona",
                "Reference": "-",
            },
        ),
    params:
        extra="-m 3 -t 5",
        tax_dir=lambda wildcards, input: os.path.dirname(input.tax_tab),
    log:
        "logs/kraken/krona_chart/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "ktImportTaxonomy {params.extra} -tax {params.tax_dir} -o {output} {input.kraken_output} 1> {log} 2>&1"
