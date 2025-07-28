"""
name: run MCScanX
date: 20234-03-22
description: Snakemake pipeline to run MCScanX
author: Jianmin Xie
dependencies:
    - seqtk = 1.4
    - blast = 2.15.0
    - MCScanX = 1.0.0
"""
import os
from os import path
import pandas as pd
shell.executable('/bin/zsh') # set a shell to run commands

configfile: "configs/MCScanX_config.yaml"

env_name = config["env_name"]
pgap_output = config["pgap_output"]

combine_faa_file = '/home/xjm/final_representative_strain/Pha/TG_genome/my_pangenome/pan_genome_output/seq/combined-aas.faa'

blast_db_output = "blast_db"
diamond_db_output = "diamond_db"
external_gene_calls_output = "external_gene_calls"
contigs_db_output = "contigs_db"
pan_genome_output = "pan_genome_output"

MCScanX_output = "output"



SAMPLES =  os.listdir(pgap_output)


rule all:
    input:
        # expand(MCScanX_output + "/gff/{sample}.gff", sample=SAMPLES)
        MCScanX_output+ "/combined.gff",
        # MCScanX_output + "/combined.faa",
        # blast_db_output + "/combined.phr",
        MCScanX_output + "/combined.blast",
        MCScanX_output + "/combined.collinearity",



rule make_gff:
    input:
        gff=pgap_output + "/{sample}/result/annot_with_genomic_fasta.gff",
    output:
        gff=MCScanX_output + "/gff/{sample}.gff",
    shell:
        """grep '\CDS\s'  {input.gff} | awk -F '[\t;=]' 'BEGIN{{OFS="\t"}} {{print $1,$10,$4,$5}}' | grep 'chromosome' | sed  's/cds-//g' > {output.gff}"""

rule combine_gff:
    input:
        expand(MCScanX_output + "/gff/{sample}.gff", sample=SAMPLES)
    output:
        gff=MCScanX_output + "/combined.gff",
        id=MCScanX_output + "/id.txt",
    shell:
        'cat {input} > {output.gff};'
        'cut -f2 {output.gff} > {output.id}'



rule combine_faa:
    input:
        id = rules.combine_gff.output.id,
        faa = combine_faa_file,
    output:
        faa=MCScanX_output + "/combined.faa",
    conda: env_name
    shell:
        'seqtk subseq {input.faa} {input.id} > {output}'




rule make_blast_db:
    input:
        rules.combine_faa.output.faa
    output:
        blast_db_output + "/combined.phr"
    params:
        db_prifix = blast_db_output + "/combined"
    conda:
        env_name
    log:
        "logs/make_blast_db.log"
    shell:
        "makeblastdb -in {input} -dbtype prot -out {params.db_prifix} &> {log}"

rule run_blastp:
    input:
        phr = rules.make_blast_db.output,
        faa = rules.combine_faa.output.faa
    output:
        MCScanX_output + "/combined.blast"
    params:
        db_prifix = blast_db_output + "/combined",
        max_target_seqs = 500, # Ensure this value > number of genomes
    conda:
        env_name
    threads: workflow.cores
    log:
        "logs/blastp.log"
    shell:
        "blastp -query {input.faa} -out {output} -db {params.db_prifix} -outfmt 6 -evalue 1e-5 -max_target_seqs {params.max_target_seqs} -num_threads {threads} &> {log}"

rule run_MCScanX:
    input:
        MCScanX_output + "/combined.blast"
    output:
        MCScanX_output + "/combined.collinearity"
    params:
        output_dir = MCScanX_output
    log:
        "../logs/MCScanX.log"
    shell:
        "cd {params.output_dir};MCScanX combined -m 5 &> {log}" # -m The maximum number of gene intervals allowed in a block (the maximum number of genes that can be between two connected genes)
