from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


pepfile: config["pepfile"]


validate(pep.sample_table, "../schemas/samples.schema.yaml")


def get_sample_names():
    return list(pep.sample_table["sample_name"].values)


def get_one_fastq_file(wildcards, read_pair="fq1"):
    return pep.sample_table.loc[wildcards.sample][[read_pair]]


def get_fastq_paths(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def get_constraints():
    return {
        "sample": "|".join(get_sample_names()),
    }


def get_outputs():
    sample_names = get_sample_names()
    outputs = {
        "fastqc_report": expand(
            "results/reads/deduplicated/fastqc/{sample}_R{orientation}.html",
            sample=sample_names,
            orientation=[1, 2],
        ),
        "kronas": expand("results/kraken/kronas/{sample}.html", sample=sample_names),
        "quast": expand("results/quast/{sample}/report.html", sample=sample_names),
    }
    return outputs


def infer_read_path(wildcards):
    if wildcards.step != "original":
        return "results/reads/{step}/{sample}_{orientation}.fastq.gz"
    else:
        if wildcards.orientation == "R1":
            return get_one_fastq_file(wildcards, read_pair="fq1")[0]
        elif wildcards.orientation == "R2":
            return get_one_fastq_file(wildcards, read_pair="fq2")[0]


#### COMMON STUFF #################################################################


def get_cutadapt_extra() -> list[str]:
    args_lst = []
    if config["reads__trimming"].get("keep_trimmed_only", False):
        args_lst.append("--discard-untrimmed")
    if "shorten_to_length" in config["reads__trimming"]:
        args_lst.append(f"--length {config['reads__trimming']['shorten_to_length']}")
    if "cut_from_start" in config["reads__trimming"]:
        args_lst.append(f"--cut {config['reads__trimming']['cut_from_start']}")
    if "cut_from_end" in config["reads__trimming"]:
        args_lst.append(f"--cut -{config['reads__trimming']['cut_from_end']}")
    if "max_n_bases" in config["reads__trimming"]:
        args_lst.append(f"--max-n {config['reads__trimming']['max_n_bases']}")
    if "max_expected_errors" in config["reads__trimming"]:
        args_lst.append(f"--max-expected-errors {config['reads__trimming']['max_expected_errors']}")
    if param_value := config["reads__trimming"].get("anywhere_adapter", ""):
        args_lst.append(f"--anywhere file:{param_value}")
    if param_value := config["reads__trimming"].get("front_adapter", ""):
        args_lst.append(f"--front file:{param_value}")
    if param_value := config["reads__trimming"].get("regular_adapter", ""):
        args_lst.append(f"--adapter file:{param_value}")
    return args_lst


def parse_paired_cutadapt_param(pe_config, param1, param2, arg_name) -> str:
    if param1 in pe_config:
        if param2 in pe_config:
            return f"{arg_name} {pe_config[param1]}:{pe_config[param2]}"
        else:
            return f"{arg_name} {pe_config[param1]}:"
    elif param2 in pe_config:
        return f"{arg_name} :{pe_config[param2]}"
    return ""


def parse_cutadapt_comma_param(config, param1, param2, arg_name) -> str:
    if param1 in config:
        if param2 in config:
            return f"{arg_name} {config[param2]},{config[param1]}"
        else:
            return f"{arg_name} {config[param1]}"
    elif param2 in config:
        return f"{arg_name} {config[param2]},0"
    return ""


def get_cutadapt_extra_pe() -> str:
    args_lst = get_cutadapt_extra()

    cutadapt_config = config["reads__trimming"]
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "max_length_r1", "max_length_r2", "--maximum-length"):
        args_lst.append(parsed_arg)
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "min_length_r1", "min_length_r2", "--minimum-length"):
        args_lst.append(parsed_arg)
    if qual_cut_arg_r1 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "--quality-cutoff"
    ):
        args_lst.append(qual_cut_arg_r1)
    if qual_cut_arg_r2 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "-Q"
    ):
        args_lst.append(qual_cut_arg_r2)
    return " ".join(args_lst)


def get_kraken_decontamination_params():
    extra = []
    if config["reads__decontamination"]["exclude_children"]:
        extra.append("--include-children")
    if config["reads__decontamination"]["exclude_ancestors"]:
        extra.append("--include-parents")
    return " ".join(extra)


def get_quast_params():
    mincontig_param = "--min-contig {val}".format(val=config["quast__params"]["min_contig_length"])
    if config["quast__params"]["extra"]:
        return f'{mincontig_param} {config["quast__params"]["extra"]}'
    return mincontig_param


def get_spades_mode():
    return "" if config["spades__params"]["mode"] == "standard" else f'--{config["spades__params"]["mode"]}'


### RESOURCES


def get_mem_mb_for_trimming(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["trimming_mem_mb"] * attempt)


def get_mem_mb_for_spades(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["spades_mem_mb"] * attempt)


def get_mem_mb_for_fastqc(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["fastqc_mem_mb"] * attempt)
