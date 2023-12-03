rule spades__assemble_reads_into_contigs:
    input:
        r1="results/reads/deduplicated/{sample}_R1.fastq.gz",
        r2="results/reads/deduplicated/{sample}_R2.fastq.gz",
    output:
        fastg="results/assembly/{sample}/assembly_graph.fastg",
        gfa="results/assembly/{sample}/assembly_graph_with_scaffolds.gfa",
        contigs="results/assembly/{sample}/contigs.fasta",
        scaffolds="results/assembly/{sample}/scaffolds.fasta",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.contigs),
        mode=get_spades_mode(),
        careful="--careful" if config["spades__params"]["careful"] else "",
    threads: min(config["threads"]["spades"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_spades,
    log:
        "logs/spades/{sample}.log",
    conda:
        "../envs/spades.yaml"
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o {params.outdir} --threads {threads} {params.mode} {params.careful} > {log} 2>&1"


rule quast__evaluate_assembly:
    input:
        fasta="results/assembly/{sample}/contigs.fasta",
    output:
        pdf=report(
            "results/quast/{sample}/report.pdf",
            category="{sample}",
            labels={"Type": "QUAST"},
        ),
        basic_reports=multiext("results/quast/{sample}/report.", "html", "tex", "txt", "tsv"),
        transposed=multiext("results/quast/{sample}/transposed_report.", "tex", "txt", "tsv"),
        basic_stats=directory("results/quast/{sample}/basic_stats/"),
        icarus="results/quast/{sample}/icarus.html",
        viewer="results/quast/{sample}/icarus_viewers/contig_size_viewer.html",
        log="results/quast/{sample}/quast.log",
    log:
        "logs/quast/{sample}.log",
    threads: min(config["threads"]["quast"], config["max_threads"])
    params:
        extra=get_quast_params(),
    wrapper:
        "v2.13.0/bio/quast"
