"""
name: extract 16S
date: 2023-11-17
description: Snakemake pipeline to extract 16S from genome
author: Jianmin Xie
dependencies:
    - barrnap = 0.9
    - biopython = 1.81
    - blast = 2.15.0
    - pandas = 2.2.3
"""
import os
from os import path
import pandas as pd

shell.executable('/bin/zsh') # set a shell to run commands

# base_dir = "/home/xjm/cyr/Germany_novel_strain/Verrucomicrobiota_genome_all"
genome_dir = config["fasta"]
suffix = config["suffix"]
output_16S_dir = "barrnap_output"
blast_db_output = "blast_db"
env_name = "barrnap"



SAMPLES = [".".join(i.split(".")[:-1]) for i in os.listdir(genome_dir) if not i.startswith(".")]
# print(len(SAMPLES))

rule all:
    input:
        "16S_identity_matrix.tsv"



rule barrnap:
    input:
        genome_dir + "/{sample}." + suffix
    output:
        fasta=output_16S_dir + "/{sample}.fasta",
        fai=temp(genome_dir + "/{sample}."+ suffix +".fai")
    log:
        "log/{sample}.log"
    conda:
        env_name
    shell:
        "barrnap  --outseq {output.fasta} {input} &> {log}"

rule extract:
    input:
        expand(output_16S_dir + "/{sample}.fasta",sample=SAMPLES)
    output:
        "extract_16S.fasta"
    conda:
        env_name
    log:
        "log/extract.log"
    script:
        "scripts/python/extract_16S_rRNA_gene.py"


rule make_blast_db:
    input:
        rules.extract.output
    output:
        directory(blast_db_output)
    log:
        "log/make_blast_db.log"
    conda:
        env_name
    shell:
        "mkdir {output};"
        "makeblastdb -in {input}  -dbtype nucl -out {output}/16S &> {log}"

rule blast:
    input:
        db=rules.make_blast_db.output,
        query=rules.extract.output
    output:
        "blast_out.tsv"
    params:
        max_target_count = len(SAMPLES)
    conda:
        env_name
    shell:
        "blastn -db {input.db}/16S -query {input.query} "
        "-outfmt '6 qacc sacc pident qcovhsp qlen slen length mismatch gapopen evalue bitscore' "
        "-out {output} -evalue 1e-5 -max_target_seqs {params.max_target_count}"

rule table_to_matrix:
    input:
        rules.blast.output
    output:
        "16S_identity_matrix.tsv"
    run:
        df = pd.read_table(input[0], header=None)
        df.columns = ["qacc", "sacc", "pident", "qcovhsp", "qlen", "slen", "length", "mismatch", "gapopen", "evalue", "bitscore"]
        df_matrix = df.pivot(index='qacc',columns='sacc',values='pident')
        df_matrix.to_csv(output[0], sep="\t")


