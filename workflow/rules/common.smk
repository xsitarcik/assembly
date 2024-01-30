from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


### Layer for adapting other workflows  ###############################################################################


def get_fastq_for_assembly(wildcards):
    reads = reads_workflow.get_final_fastq_for_sample(wildcards.sample)
    return {
        "r1": reads[0],
        "r2": reads[1],
    }


def get_sample_names():
    return reads_workflow.get_sample_names()


### Data input handling independent of wildcards ######################################################################


def parse_blast_tag(path: str):
    return os.path.basename(os.path.dirname(os.path.realpath(path)))


def validate_blast_tag(tag: str):
    VALID_TAGS = [
        "18S_fungal_sequences",
        "Betacoronavirus",
        "ref_viroids_rep_genomes",
        "ref_viruses_rep_genomes",
        "refseq_select_rna",
        "refseq_select_prot",
        "refseq_protein",
        "refseq_rna",
        "16S_ribosomal_RNA",
        "ITS_RefSeq_Fungi",
        "28S_fungal_sequences",
        "ITS_eukaryote_sequences",
        "LSU_eukaryote_rRNA",
        "LSU_prokaryote_rRNA",
        "SSU_eukaryote_rRNA",
        "env_nt",
        "env_nr",
        "human_genome",
        "landmark",
        "mito",
        "mouse_genome",
        "nr",
        "nt_euk",
        "nt",
        "nt_others",
        "swissprot",
        "tsa_nr",
        "tsa_nt",
        "taxdb",
        "nt_prok",
        "nt_viruses",
        "pataa",
        "patnt",
        "pdbaa",
        "pdbnt",
        "ref_euk_rep_genomes",
        "ref_prok_rep_genomes",
    ]
    if tag not in VALID_TAGS:
        raise ValueError(f"{tag=} was inferred as Blast DB tag, which is not valid. {VALID_TAGS=}")


def get_blast_ref_tag():
    return parse_blast_tag(config["assembly__classification__blast"]["db_dir"])


validate_blast_tag(get_blast_ref_tag())


### Global rule-set stuff #############################################################################################


BLAST_HEADER = "qseqid sacc staxid sscinames scomnames stitle pident evalue length mismatch gapopen qstart qend sstart send qlen slen"


def get_blast_binary():
    db_type = config["assembly__classification__blast"]["query_vs_db"]
    blast_binaries = {
        "nucleotide-nucleotide": "blastn",
        "nucleotide-protein": "blastx",
        "protein-nucleotide": "tblastn",
        "protein-protein": "blastp",
    }
    return blast_binaries[db_type]


def get_max_number_of_hits():
    return config["assembly__classification__blast"]["max_number_of_hits"]


def get_blast_db():
    blast_type = config["assembly__classification__blast"]["query_vs_db"].split("-")[1][0]
    reference_tag = get_blast_ref_tag()
    return multiext(
        os.path.join(config["assembly__classification__blast"]["db_dir"], f"{reference_tag}.{blast_type}"),
        "db",
        "ot",
        "tf",
        "to",
    )


def get_outputs():
    sample_names = get_sample_names()

    outputs = {}
    assembly_tools = config["assembly"]["assembly"]
    for tool in assembly_tools:
        outputs[tool] = expand(f"results/assembly/{{sample}}/{tool}/contigs.fasta", sample=sample_names)

    if "quast" in config["assembly"]["report"]:
        outputs["quast"] = expand(
            "results/assembly/{sample}/{assembly_tool}/QUAST/report.pdf",
            sample=sample_names,
            assembly_tool=assembly_tools,
        )
    if "bandage" in config["assembly"]["report"]:
        outputs["bandage"] = expand(
            "results/assembly/{sample}/{assembly_tool}/bandage/bandage.{ext}",
            sample=sample_names,
            assembly_tool=assembly_tools,
            ext=["info", "svg"],
        )

    if "blast" in config["assembly"]["classification"]:
        outputs["blast"] = expand(
            "results/classification/{sample}/{assembly_tool}/blast/summary.html",
            sample=sample_names,
            assembly_tool=assembly_tools,
        )

    return outputs


### Contract for other workflows ######################################################################################


### Parameter parsing from config #####################################################################################


def get_quast_params():
    mincontig_param = "--min-contig {val}".format(val=config["assembly__report__quast"]["min_contig_length"])
    if config["assembly__report__quast"]["extra"]:
        return f'{mincontig_param} {config["assembly__report__quast"]["extra"]}'
    return mincontig_param


def get_spades_params():
    mode = (
        ""
        if config["assembly__assembly__spades"]["mode"] == "standard"
        else f'--{config["assembly__assembly__spades"]["mode"]}'
    )
    careful = "--careful" if config["assembly__assembly__spades"]["careful"] else ""
    if mode and careful:
        return f"{mode} {careful}"
    return mode + careful


### Resource handling #################################################################################################


def get_mem_mb_for_assembly(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__assembly_mem_mb"] * attempt)


def get_mem_mb_for_classification(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__classification_mem_mb"] * attempt)


def get_threads_for_assembly():
    return min(config["threads"]["assembly__assembly"], config["max_threads"])


def get_threads_for_classification():
    return min(config["threads"]["assembly__classification"], config["max_threads"])


def get_threads_for_report():
    return min(config["threads"]["assembly__report"], config["max_threads"])
