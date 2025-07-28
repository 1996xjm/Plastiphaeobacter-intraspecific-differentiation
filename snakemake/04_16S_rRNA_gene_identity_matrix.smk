"""
name: 16S rRNA gene identity matrix
date: 2023-11-17
description: Snakemake pipeline to calculate 16S rRNA gene identity matrix
author: Jianmin Xie
dependencies:
    - blast = 2.15.0
    - pandas = 2.2.1
rules:
    - barrnap
"""
import os
from os import path
import pandas as pd

shell.executable('/bin/zsh') # set a shell to run commands

fasta = config["fasta"]
blast_db_output = "blast_db"
env_name = "barrnap"
max_target_count = 47



rule all:
    input:
        "16S_identity_matrix.tsv"



rule make_blast_db:
    input:
        fasta
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
        query=fasta
    output:
        "blast_out.tsv"
    params:
        max_target_count = max_target_count
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

        df_matrix = df[df["qcovhsp"]>50].pivot(index='qacc',columns='sacc',values='pident')
        df_matrix.to_csv(output[0], sep="\t")


