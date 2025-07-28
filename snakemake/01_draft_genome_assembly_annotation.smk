"""
name: Second generation sequencing assembly
date: 2023-08-27
description: Snakemake pipeline to assembly second generation sequencing reads.
author: Jianmin Xie
dependencies:
    - fastp = 0.23.4
    - spades = 3.15.5
    - biopython = 1.81
    - checkm-genome = 1.2.2
    - busco = 5.5.0
    - quast = 5.2.0
    - bakta = 1.9.2
    - eggnog-mapper = 2.1.12

"""
import os
configfile: "configs/SG_assembly_config.yaml"

sample_path = config["samples_path"]
bakta_db = config["bakta_db"]
env_name = config["env_name"]


FASTP_OUT_PATH = "fastp_out"
SPADES_OUT_PATH = "spades_out"
FILTERED_GENOME_OUT_PATH = "genome_upper_L1000C10"
BUSCO_OUT_PATH = "busco_output"
CHECKM_OUT_PATH = "checkM_output"
EGGNOG_OUT_PATH = "eggnog_out"
bakta_OUT_PATH = "bakta_output"
bakta_env_name = "bakta"


SAMPLES = os.listdir(sample_path)



module fastp_m:
    snakefile: "modules/fastp_module.smk"
    config: config["fastp_config"]

quality_assessment_config = {
    "genome_dir": FILTERED_GENOME_OUT_PATH,
    "env_name": env_name,
    "samples": SAMPLES,
    "busco_database": config["busco_database"],
    "checkM_output": CHECKM_OUT_PATH,
    "busco_output": BUSCO_OUT_PATH
}

module quality_assessment_m:
    snakefile: "modules/quality_assessment_module.smk"
    config: quality_assessment_config



use rule fastp from fastp_m as my_fastp

use rule * from quality_assessment_m as my_*


rule all:
    input:
        BUSCO_OUT_PATH + "/batch_summary.tsv",
        CHECKM_OUT_PATH+ "/checkm_result.tsv",
        expand(bakta_OUT_PATH + "/{sample}/{sample}.faa", sample=SAMPLES),
        # expand(EGGNOG_OUT_PATH + "/{sample}/{sample}.emapper.annotations", sample=SAMPLES),
        # expand(FILTERED_GENOME_OUT_PATH + "/{sample}.fasta", sample=SAMPLES),
    default_target: True

rule spades:
    input:
        r1 = rules.my_fastp.output.r1,
        r2 = rules.my_fastp.output.r2,
    threads: workflow.cores * 0.9
    resources:
        mem_gb = 96
    output:
        ensure(SPADES_OUT_PATH + "/{sample}/contigs.fasta", non_empty=True)
    params:
        out = SPADES_OUT_PATH + "/{sample}"
    benchmark:
        "benchmarks/spades/{sample}_benchmark.txt"
    conda:
        env_name
    log:
        "logs/spades/{sample}_spades.log"
    shell:
        "spades.py --isolate -1 {input.r1} -2 {input.r2} -t {threads} -m {resources.mem_gb} -o {params.out} &> {log}"



rule filter_contig:
    input:
        rules.spades.output
    output:
        FILTERED_GENOME_OUT_PATH + "/{sample}.fasta"
    params:
        thLength = 1000,  # Length threshold
        thCoverage = 10,  # Coverage threshold
    log:
        "logs/filter/{sample}_filter_contig.log"
    conda:
        env_name
    script:
        "scripts/python/filterLengthAndCoverage.py"


rule bakta:
    input:
        genome = rules.filter_contig.output,
    output:
        bakta_OUT_PATH + "/{sample}/{sample}.faa"
    params:
        out_dir = bakta_OUT_PATH + "/{sample}",
        db = bakta_db,
    log:
        "logs/bakta/{sample}_bakta.log"
    threads: workflow.cores/2
    conda:
        bakta_env_name
    shell:
        "bakta -f "
        "--db {params.db} "
        "--verbose "
        "--output {params.out_dir} "
        "--prefix {wildcards.sample} "
        "--locus-tag {wildcards.sample} "
        "--threads {threads} {input.genome} "
        "&> {log}"

rule eggnog:
    input:
        rules.bakta.output
    output:
        EGGNOG_OUT_PATH + "/{sample}/{sample}.emapper.annotations"
    threads: workflow.cores * 0.6
    params:
        prefix = EGGNOG_OUT_PATH + "/{sample}/{sample}",
        db_dir = config["eggnog_db"]
    log:
        "logs/eggnog/{sample}_eggnog.log"
    conda:
        env_name
    shell:
        "emapper.py --override --cpu {threads} -i {input} -o {params.prefix} --data_dir {params.db_dir} &> {log}"
