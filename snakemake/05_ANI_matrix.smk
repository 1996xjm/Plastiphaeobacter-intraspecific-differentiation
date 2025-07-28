"""
name: run fastANI
date: 2023-04-07
description: Snakemake pipeline to run fastANI
author: Jianmin Xie
dependencies:
    - fastANI = 1.33
    - pandas = 2.2.3
"""
import os
from os import path
from io import StringIO
import pandas as pd
import numpy as np
shell.executable('/bin/zsh') # set a shell to run commands

genome_dir = config["fasta"]
fastANI_output = "fastANI_output"

genome_abs_dir = path.abspath(genome_dir)

SAMPLE = ["Pha"]



rule all:
    input:
        "ANI_matrix.tsv"



rule make_id_text:
    input:
        genome_dir
    output:
        ql = fastANI_output + "/ql.txt",
        rl = fastANI_output + "/rl.txt"
    run:
        fd_ql = open(output.ql, "w")
        fd_rl = open(output.rl, "w")
        for sample in os.listdir(input[0]):
            fd_ql.write(path.abspath(path.join(input[0],sample))+"\n")
            fd_rl.write(path.abspath(path.join(input[0],sample))+"\n")
        fd_ql.close()
        fd_rl.close()

rule fastANI:
    input:
        ql = rules.make_id_text.output.ql,
        rl = rules.make_id_text.output.rl
    output:
        table = fastANI_output + "/ANI.out",
        matrix = fastANI_output + "/ANI.out.matrix"
    params:
        rep = genome_abs_dir
    conda:
        "fastani"
    threads: 32
    shell:
        "fastANI --ql {input.ql} --rl {input.rl} --matrix -o {output.table} -t {threads};"
        "sed -i 's#{params.rep}/##g' {output.matrix};"
        "sed -i 's#.fasta##g' {output.matrix};"


rule table_to_matrix:
    input:
        rules.fastANI.output.matrix
    output:
        "ANI_matrix.tsv"
    run:
        with open(input[0],"r") as f:
            string_io = StringIO()
            line_num = 0
            for line in f:
                if line_num == 0:
                    sample_count = int(line.strip())
                    string_io.write("\t".join(["sample_id"]+[str(i) for i in range(sample_count)])+"\n")
                else:
                    string_io.write(line)
                line_num += 1
            string_io.seek(0)
            df = pd.read_table(string_io, index_col=0)
            # print(df)
            # Copy the values of the lower triangular matrix to the upper triangular
            whole_matrix = np.tril(df.values) + np.triu(df.T.values)
            matrix_df = pd.DataFrame(whole_matrix,columns=df.index,index=df.index)
            
            for i in range(len(matrix_df)):
                if pd.isna(matrix_df.iat[i, i]):
                    matrix_df.iat[i, i] = 100

            
            matrix_df = matrix_df.fillna(0)
            matrix_df.to_csv(output[0], sep="\t")

