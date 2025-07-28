"""
name: The pipeline to build MSA genome tree
date: 2023-04-24
description: Snakemake pipeline to do Third generation sequence assembly
author: Jianmin Xie
dependencies:
    - prodigal = 2.6.3
    - biopython = 1.81
    - orthofinder = 2.5.5
    - mafft = 7.520
    - seqtk = 1.4
    - iqtree = 2.2.3
"""

import os

configfile: "configs/MSA_genome_tree_config.yaml"

genomes_dir = config["genomes_dir"]
env_name = config["env_name"]
SAMPLES = [".".join(i.split(".")[:-1]) for i in os.listdir(genomes_dir)]

paired_count = len(SAMPLES) ** 2
ulimit_count = 10000 if paired_count < 10000 else paired_count







prodigal_OUT_PATH = "prodigal_output"
seqtk_OUT_PATH = "seqtk_output"
orthofinder_OUT_PATH = "orthofinder_output"
scos_dir = "Single_Copy_Orthologue_Sequences"
mafft_dir = "mafft_output"
iqtree_dir = "iqtree_output"

shell.executable('/bin/zsh') # set a shell to run commands

def get_OG_names(wildcards):
    # note 1: ck_output is the same as OUTDIR, but this requires
    # the checkpoint to complete before we can figure out what it is!

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name. Here we use make_some_files
    filter_output = checkpoints.filterSingleCopyGene.get(**wildcards).output[0]
    OGs, = glob_wildcards(os.path.join(filter_output, "{OG}.fa"))
    return expand(os.path.join(mafft_dir, "{OG}.faa"), OG=OGs)

rule all:
    input:
        scos_dir,
        iqtree_dir + "/MSA_genome_tree.treefile"
    default_target: True


# mapping
rule prodigal:
    input:
        genomes_dir + "/{sample}.fasta"
    output:
        temp(prodigal_OUT_PATH + "/{sample}.faa")

    log:
        "logs/prodigal/{sample}_flye.log"

    conda:
        env_name
    shell:
        "prodigal -i {input} -a {output} &> {log}"


rule seqtk:
    input:
        rules.prodigal.output
    output:
        seqtk_OUT_PATH + "/{sample}.faa"

    params:
        prefix = "{sample}_"
    conda:
        env_name
    shell:
        "seqtk rename {input} {params.prefix} > {output}"



rule orthofinder:
    input:
        expand(seqtk_OUT_PATH + "/{sample}.faa",sample=SAMPLES)
    output:
        ensure(orthofinder_OUT_PATH + "/Results_defaut/Comparative_Genomics_Statistics/Statistics_Overall.tsv", non_empty=True)

    params:
        faa_dir = seqtk_OUT_PATH,
        out_dir = orthofinder_OUT_PATH,
        ulimit_count = ulimit_count
    log:
        "logs/orthofinder.log"
    conda:
        env_name
    shell:
        "rm -rf {params.out_dir};"
        "ulimit -n {params.ulimit_count};"
        "orthofinder -f {params.faa_dir}  -o {params.out_dir} -n defaut &> {log}"



checkpoint filterSingleCopyGene:
    input:
        rules.orthofinder.output
    output:
        directory(scos_dir)
    params:
        single_copy_gene_dir = orthofinder_OUT_PATH + "/Results_defaut/Single_Copy_Orthologue_Sequences"
    log:
        "logs/filterSingleCopyGene.log"
    conda:
        env_name
    script:
        "scripts/python/filterSingleCopyGene.py"




rule mafft:
    input:
        scos_dir + "/{OG}.fa"
    output:
        mafft_dir + "/{OG}.faa"
    log:
        "logs/mafft/{OG}.log"
    conda:
        env_name
    shell:
        "(mafft --auto {input} > {output}) &> {log}"





rule concatenate_seq:
    input:
        get_OG_names
    output:
        iqtree_dir + "/concatenated_seq.fasta"
    conda:
        env_name
    script:
        "scripts/python/concatenate_seq.py"




rule iqtree:
    input:
        rules.concatenate_seq.output
    output:
        iqtree_dir + "/MSA_genome_tree.treefile"
    params:
        outgroup = config["outgroup"] if "outgroup" in config else "\"\"",
        output_prefix = iqtree_dir + "/MSA_genome_tree"
    log:
        "logs/iqtree.log"
    conda:
        env_name
    threads: workflow.cores
    shell:
        """if [[ {params.outgroup} == "" ]] {{
    iqtree2 -s {input} -m MFP  -B 1000  -redo --prefix {params.output_prefix}  -nt {threads} &> {log}
}} else {{
    iqtree2 -s {input} -m MFP  -B 1000  -redo --prefix {params.output_prefix}  -nt {threads} -o {params.outgroup} &> {log}
}}"""













