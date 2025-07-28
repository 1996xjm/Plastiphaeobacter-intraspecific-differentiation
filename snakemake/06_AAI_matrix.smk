"""
name: run EzAAI
date: 2023-02-27
description: Snakemake pipeline to run EzAAI
author: Jianmin Xie
dependencies:
    - EzAAI = 1.2.2
"""
import os
from os import path
import pandas as pd
shell.executable('/bin/zsh') # set a shell to run commands

fasta_dir = config["fasta"]
env_name = config["env"]
output_DB_dir = "DBs"
output_file = "AAI.tsv"
output_matrix_file = "AAI_matrix.tsv"


SAMPLE = [".".join(i.split(".")[:-1]) for i in os.listdir(fasta_dir) if not i.startswith(".")]

onsuccess:
    df = pd.read_table(output_file)
    matrix = df["AAI"].groupby(by=[df["Label 1"], df["Label 2"]]).sum().unstack()
    matrix.to_csv(output_matrix_file,sep="\t")
rule all:
    input:
        output_file


rule make_db:
    params:
        label="{sample}"
    input:
        fasta_dir + "/{sample}.fasta"
    output:
        output_DB_dir + "/{sample}.msd"
    log:
        "logs/make_db/{sample}.log"
    conda:
        env_name
    shell:
        # "EzAAI convert -i {input} -s prot -o {output} -l {params.label} &> {log}" # Protein sequences
        'EzAAI extract -i "{input}" -o "{output}" -l "{params.label}" &> "{log}"' # Nucleotide sequences


rule EzAAI:
    input:
        expand(output_DB_dir + "/{sample}.msd",sample=SAMPLE)
    output:
        output_file
    params:
        db_dir = output_DB_dir
    log:
        "logs/EzAAI.log"
    threads: 128
    conda:
        env_name
    shell:
        "EzAAI calculate -i {params.db_dir} -j {params.db_dir} -o {output} -t {threads} &> {log}"