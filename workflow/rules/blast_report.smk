rule custom__split_contigs:
    input:
        contigs="results/assembly/{sample}/contigs.fasta",
    output:
        fasta_dir=directory("results/summary_report/{sample}/annotation/sequences/"),
        seqinfo="results/summary_report/{sample}/annotation/attributes/seqinfo.fa",
    params:
        fasta=lambda wildcards, output: os.path.join(output.fasta_dir, "all_contigs.fa"),
    log:
        "logs/custom/split_contigs/{sample}.log",
    conda:
        "../envs/blast_report.yaml"
    script:
        "../scripts/summary_contigs.py"


rule custom__summary_blast:
    input:
        tsv="results/blast/{sample}/{reference_tag}.tsv",
    output:
        tsv="results/summary_report/{sample}/annotation/attributes/blast/{reference_tag}.blast.tsv",
    params:
        min_query_coverage=config["blast__report"]["min_query_coverage"],
        max_target_seqs=config["blast__report"]["max_target_seqs"],
    log:
        "logs/custom/summary_blast/{sample}.log",
    conda:
        "../envs/blast_report.yaml"
    script:
        "../scripts/summary_blast.py"


rule custom__summary_html:
    input:
        template="../resources/abundances.html",
        fasta_dir="results/summary_report/{sample}/annotation/sequences/",
        attrs_seqinfo="results/summary_report/{sample}/annotation/attributes/seqinfo.fa",
        attrs_blast=get_all_blast_results,
    output:
        html="results/summary_report/{sample}/summary.html",
        tsv="results/summary_report/{sample}/summary.html.tsv",
    params:
        fasta=lambda wildcards, input: os.path.join(input.fasta_dir, "all_contigs.fa"),
        max_query_seqs=config["blast__report"]["max_query_seqs"],
        seqs_per_page=config["blast__report"]["seqs_per_page"],
        sort_by="Sequence",
        sort_how="asc",
        fasta_dir="sequences",
        html_columns=["Sequence", "Length", "Compress ratio", "Homologue link"],
    log:
        "logs/custom/summary_html/{sample}.log",
    conda:
        "../envs/blast_report.yaml"
    script:
        "../scripts/summary.py"
