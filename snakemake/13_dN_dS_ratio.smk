"""
name: Use CODEML to estimate the synonymous and nonsynonymous rates (dS and dN) and to detect positive Darwinian selection
date: 2024-05-24
description: Snakemake pipeline for The CODEML program in the PAML package has been widely used to analyze protein-coding gene sequences to estimate the synonymous and nonsynonymous rates (dS and dN) and to detect positive Darwinian selection driving protein evolution.
reference: https://doi.org/10.1093/molbev/msad041, https://github.com/abacus-gene/paml-tutorial/tree/main/positive-selection
author: Jianmin Xie
dependencies:
    - pal2nal.pl http://www.bork.embl.de/pal2nal/
    - mafft = 7.520
    - raxml-ng = 1.2.0
    - codeml in paml = 4.10.7
    - pandas = 2.2.3

"""

import os
from os import path
import pandas as pd
from Bio import SeqIO, AlignIO
from io import StringIO
import re



shell.executable('/bin/zsh') # set a shell to run commands


single_copy_gene_fnas_dir = "/home/xjm/final_representative_strain/Pha/TG_genome/my_pangenome/pan_genome_output/seq/single_copy_fnas"
single_copy_gene_faas_dir = "/home/xjm/final_representative_strain/Pha/TG_genome/my_pangenome/pan_genome_output/seq/single_copy_aas"
SAMPLES = [".".join(i.split(".")[:-1]) for i in os.listdir(single_copy_gene_fnas_dir) if i.endswith("fna")]
single_copy_gene_info_file = "/home/xjm/final_representative_strain/Pha/TG_genome/my_pangenome/pan_genome_output/table/single_copy_genes_annotation/PT42_8_single_copy_genes_annotation.tsv"
alignments_output = "alignments_output"
good_SNV_density_alignments_output = "good_SNV_density_alignments_output"
raxml_output = "raxml_output"
codeml_output = "codeml_output"
summary_output = "summary_output"

env_name = "orthofinder"


def get_GC_names(wildcards):
    # note 1: ck_output is the same as OUTDIR, but this requires
    # the checkpoint to complete before we can figure out what it is!

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name. Here we use make_some_files
    filter_output = checkpoints.filter_GC_by_SNV_density.get(**wildcards).output[0]
    GCs, = glob_wildcards(os.path.join(filter_output, "{GC}_aln.fna"))
    return expand(codeml_output + "/branch_free_ratio_models/{GC}/out_branch.txt", GC=GCs)

rule all:
    input:
        summary_output+ "/branch_free_ratio_models/omega_matrix.tsv",
        summary_output+ "/density/SNV_density_info.tsv",
        expand(codeml_output + "/branch_free_ratio_models/{sample}/out_branch.txt", sample=SAMPLES),
        # expand(codeml_output + "/homogenous_model/{sample}/out_M0.txt", sample=SAMPLES),
        # expand(codeml_output + "/site_models/{sample}/LRT_site_models.png", sample=SAMPLES),
        # summary_output + "/branch_free_ratio_models/omega_matrix.tsv",






rule mafft_faa:
    input:
        single_copy_gene_faas_dir + "/{sample}.faa"
    output:
        alignments_output + "/{sample}/{sample}_aln.faa"
    conda:
        env_name
    shell:
        "mafft --quiet --auto {input} > {output}"



rule pal2nal:
    input:
        fna = single_copy_gene_fnas_dir + "/{sample}.fna",
        aln_faa = rules.mafft_faa.output
    output:
        fasta = ensure(alignments_output + "/{sample}/{sample}_aln.fna", non_empty=True), # 如果下面的命令转换失败，输出的文件就是空白
        paml = ensure(alignments_output + "/{sample}/{sample}_aln.paml", non_empty=True),
    shell:
        "pal2nal.pl {input.aln_faa} {input.fna} -output fasta -codontable 11 -nogap > {output.fasta};"
        "pal2nal.pl {input.aln_faa} {input.fna} -output paml -codontable 11 -nogap > {output.paml}"


rule get_SNV_SAAV_density:
    input:
        fna_list = lambda x:[alignments_output + f"/{sample}/{sample}_aln.fna" for sample in SAMPLES if path.exists(alignments_output + f"/{sample}/{sample}_aln.fna")],
        faa_list = lambda x:[alignments_output + f"/{sample}/{sample}_aln.faa" for sample in SAMPLES if path.exists(alignments_output + f"/{sample}/{sample}_aln.fna")],
        single_copy_gene_info_file= single_copy_gene_info_file
    output:
        SNV_density = summary_output+ "/density/SNV_density_info.tsv",
        SAAV_density = summary_output+ "/density/SAAV_density_info.tsv",
        merge_density = summary_output+ "/density/merge_density.tsv",
    log:
        "logs/get_SNV_SAAV_density.log"
    script:
        "scripts/python/get_SNV_SAAV_density_of_seq_alignment.py"








rule raxml:
    input:
        rules.pal2nal.output.fasta
    output:
        raxml_output + "/{sample}/{sample}.raxml.bestTree"
    params:
        prefix = raxml_output + "/{sample}/{sample}"
    log:
        "logs/raxml/{sample}.log"
    conda:
        env_name
    threads: 16
    shell:
        """raxml-ng --all --msa {input} --model GTR+G --prefix {params.prefix} --threads {threads} --seed 2 --redo &> {log}
speciesNum=`grep '>' {input}|wc -l`
sed -i "1i\   $speciesNum  1" {output}"""

"""
Here, we use the mutation-selection model with observed codon frequencies used as estimates (CodonFreq=7, estFreq=0) (Yang and Nielsen 2008). 
This model explicitly accounts for the mutational bias and selection affecting codon usage, and is preferable over the other models concerning codon usage (Yang and Nielsen 2008). 
"""





rule make_ctl_file_branch_free_ratio_models:
    input:
        seqfile= alignments_output + "/{sample}/{sample}_aln.paml",
        treefile = raxml_output + "/{sample}/{sample}.raxml.bestTree"
    output:
        codeml_output + "/branch_free_ratio_models/{sample}/codeml-branch.ctl"
    run:
        ctl_str = f"""seqfile = {path.abspath(input["seqfile"])}           * Path to the alignment file
treefile = {path.abspath(input["treefile"])}           * Path to the tree file
outfile = out_branch.txt            * Path to the output file

noisy = 3              * How much rubbish on the screen
verbose = 1              * More or less detailed report

seqtype = 1              * Data type
ndata = 1           * Number of data sets or loci
icode = 0              * Genetic code
cleandata = 0              * Remove sites with ambiguity data?

model = 1         * 1 means one ratio for each branch (the free-ratio model)
NSsites = 0          * Models for ω varying across sites
CodonFreq = 7        * mutation-selection model
estFreq = 0        * Use observed freqs or estimate freqs by ML
clock = 0          * Clock model
fix_omega = 0         * Estimate or fix omega
omega = 0.5        * Initial or fixed omega
"""
        with open(output[0], "w") as f:
            f.write(ctl_str)


rule codeml_branch_free_ratio_models:
    input:
        rules.make_ctl_file_branch_free_ratio_models.output
    output:
        codeml_output + "/branch_free_ratio_models/{sample}/out_branch.txt"
    log:
        "logs/codeml_branch_free_ratio_models/{sample}.log"
    params:
        pwd = codeml_output + "/branch_free_ratio_models/{sample}"
    shell:
        "(cd {params.pwd};codeml codeml-branch.ctl) &> {log}"




rule summary_omega_from_branch_free_ratio_models:
    input:
        codeml_out_list = lambda x:[codeml_output + f"/branch_free_ratio_models/{i}/out_branch.txt" for i in SAMPLES if path.exists(codeml_output + f"/branch_free_ratio_models/{i}/out_branch.txt")],
        SNV_density = summary_output+ "/density/SNV_density_info.tsv",
        single_copy_gene_info_file = single_copy_gene_info_file
    output:
        omega_matrix = summary_output + "/branch_free_ratio_models/omega_matrix.tsv",
        omega_PCA = summary_output + "/branch_free_ratio_models/omega_log10_PCA.tsv",

        omega_SNV_density_greater_10percent_PCA = summary_output + "/branch_free_ratio_models/omega_SNV_density_greater_10percent_log10_PCA.tsv",
        omega_SNV_density_greater_10percent_matrix_log10 = summary_output + "/branch_free_ratio_models/omega_SNV_density_greater_10percent_log10_matrix.tsv",
        omega_SNV_density_greater_10percent_matrix = summary_output + "/branch_free_ratio_models/omega_SNV_density_greater_10percent_matrix.tsv",
    conda:
        "call_SNPs"
    log:
        "logs/summary_omega_from_branch_free_ratio_models.log"
    script:
        "scripts/python/summary_omega_from_branch_free_ratio_models.py"

