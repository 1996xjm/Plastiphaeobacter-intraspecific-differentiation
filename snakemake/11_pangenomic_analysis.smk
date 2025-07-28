"""
name: run anvi pangenomics
date: 2023-03-22
description: Snakemake pipeline to run anvi pangenomics https://merenlab.org/2016/11/08/pangenomics-v2/
author: Jianmin Xie
dependencies:
    - anvio = 8
    - blast = 2.15.0
"""
import os
from os import path
import pandas as pd
shell.executable('/bin/zsh') # set a shell to run commands

configfile: "configs/anvio_pangenome_config.yaml"

env_name = config["env_name"]
pgap_output = config["pgap_output"]
merge_annotation_output = config["merge_annotation_output"]


blast_db_output = "blast_db"
diamond_db_output = "diamond_db"
external_gene_calls_output = "external_gene_calls"
contigs_db_output = "contigs_db"
pan_genome_output = "pan_genome_output"



SAMPLES =  os.listdir(pgap_output)


rule all:
    input:
        pan_genome_output + "/table/gene_cluster_count.tsv"



rule rename_and_merge_pgap_faa:
    input:
        faa=expand(pgap_output + "/{sample}/result/annot_translated_cds.faa", sample=SAMPLES),
        fna=expand(pgap_output + "/{sample}/result/annot_cds_from_genomic.fna", sample=SAMPLES),
    output:
        faa=pan_genome_output + "/seq/combined-aas.faa",
        fna=pan_genome_output + "/seq/combined-cds.fna",
    log:
        "logs/rename_and_merge_pgap_faa.log"
    script:
        "scripts/python/rename_and_merge_pgap_faa.py"


rule make_blast_db:
    input:
        rules.rename_and_merge_pgap_faa.output.faa
    output:
        blast_db_output + "/combined-aas.phr"
    params:
        db_prifix = blast_db_output + "/combined-aas"
    conda:
        env_name
    log:
        "logs/make_blast_db.log"
    shell:
        "makeblastdb -in {input} -dbtype prot -out {params.db_prifix} &> {log}"

rule run_blastp:
    input:
        phr = rules.make_blast_db.output,
        faa = rules.rename_and_merge_pgap_faa.output.faa
    output:
        pan_genome_output + "/table/blast-search-results.tsv"
    params:
        db_prifix = blast_db_output + "/combined-aas",
        max_target_seqs = 10000, # Ensure this value > number of genomes
    conda:
        env_name
    threads: workflow.cores
    log:
        "logs/blastp.log"
    shell:
        "blastp -query {input.faa} -db {params.db_prifix} -evalue 1e-05 -outfmt '6 qseqid sseqid pident qcovhsp qlen slen length mismatch gapopen evalue bitscore' -max_target_seqs {params.max_target_seqs} -out {output} -num_threads {threads} &> {log}"


rule generate_mcl_input:
    input:
        rules.run_blastp.output
    output:
        pan_genome_output + "/table/mcl-input.tsv"
    params:
        minbit = 0.5 # minbit = BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B)), Threshold for filtering out dissimilar edges
    log:
        "logs/generate_mcl_input.log"
    script:
        "scripts/python/pangenome/generate_mcl_input.py"



rule run_mcl:
    input:
        rules.generate_mcl_input.output
    output:
        pan_genome_output + "/table/mcl-clusters.tsv"
    params:
        inflation = 10
    conda:
        env_name
    threads: 64
    log:
        "logs/run_mcl.log"
    shell:
        "mcl {input} --abc -I {params.inflation} -o {output} -te {threads} &> {log}"



rule generate_gene_cluster_output:
    input:
        mcl = rules.run_mcl.output[0],
        faa = rules.rename_and_merge_pgap_faa.output.faa,
        fna = rules.rename_and_merge_pgap_faa.output.fna,
        blastp_res = rules.run_blastp.output[0],
    output:
        gene_cluster_count = pan_genome_output + "/table/gene_cluster_count.tsv",
        gene_type_table = pan_genome_output + "/table/gene_type_table.tsv",
        single_copy_aas_dir = directory(pan_genome_output + "/seq/single_copy_aas"),
        single_copy_sample_aas_dir = directory(pan_genome_output + "/seq/single_copy_sample_aas"),
        all_cluster_aas_dir = directory(pan_genome_output + "/seq/all_cluster_aas"),
        single_copy_fnas_dir = directory(pan_genome_output + "/seq/single_copy_fnas"),
        single_copy_genes_annotation_dir = directory(pan_genome_output + "/table/single_copy_genes_annotation"),
        genome_specific_aas_dir = directory(pan_genome_output + "/seq/genome_specific_aas"),
        genome_specific_genes_annotation_dir = directory(pan_genome_output + "/table/genome_specific_genes_annotation"),
        genome_specific_genes_blastp_info_dir = directory(pan_genome_output + "/table/genome_specific_genes_blastp_info"),
    params:
        merge_annotation_output = merge_annotation_output
    log:
        "logs/generate_gene_cluster_output.log"
    script:
        "scripts/python/pangenome/generate_gene_cluster_table.py"

