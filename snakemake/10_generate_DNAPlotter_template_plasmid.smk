"""
name: generate DNAPlotter template
date: 2024-08-04
description: Snakemake pipeline to generate DNAPlotter (Artemis) template file.
reference:
author: Jianmin Xie
dependencies:
    - pandas = 2.2.3
    - biopython = 1.81


rules:
    -


"""
import os
from os import path
shell.executable('/bin/zsh')

plasmid_gbk = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter gallaeciensis/TG_genome/genome/plasmid_gbk'

mcl_cluster_file = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter gallaeciensis/TG_genome/my_pangenome/pan_genome_output/table/mcl-clusters.tsv'
ref_annotation_dir = f'/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter gallaeciensis/TG_genome/merge_annotation'

SAMPLES = sorted(os.listdir(plasmid_gbk))

plasmid_gff_output = "plasmid"
print(SAMPLES)
sample_list = []
gbk_list = []
for sample in SAMPLES:
    for i in os.listdir(path.join(plasmid_gbk,sample)):
        if i.endswith("gbk"):
            sample_list.append(sample)
            gbk_list.append(i.replace(".gbk",""))




rule all:
    input:
        expand(plasmid_gff_output + "/{sample}/{gbk}/orthologous_gene_track_template.tsv", zip, sample=sample_list, gbk=gbk_list),
        expand(plasmid_gff_output + "/{sample}/{gbk}/gff/transposase_genes.gff", zip, sample=sample_list, gbk=gbk_list),
        expand(plasmid_gff_output + "/{sample}/{gbk}/gff/recombinase_integrase_genes.gff", zip, sample=sample_list, gbk=gbk_list),
        expand(plasmid_gff_output + "/{sample}/{gbk}/DNAPlotter_template.txt", zip, sample=sample_list, gbk=gbk_list),




rule find_orthologous_gene:
    input:
        mcl_cluster_file = mcl_cluster_file,
        ref_annotation_file = ref_annotation_dir + "/{sample}_annotation.tsv",


    output:
        plasmid_gff_output + "/{sample}/{gbk}/orthologous_gene_track_template.tsv",
    params:
        all_smaples = SAMPLES,
        gff_out_dir = plasmid_gff_output + "/{sample}/{gbk}/gff",
        track_start_position = 0.8,
        location = "{gbk}"
    script:
        "scripts/python/DNAPlotter/find_orthologous_gene_to_gff.py"



rule find_interesting_gene:
    input:
        ref_annotation_file = ref_annotation_dir + "/{sample}_annotation.tsv",

    output:
        transposase_gene_gff_file = plasmid_gff_output + "/{sample}/{gbk}/gff/transposase_genes.gff",
        recombinase_integrase_gene_gff_file = plasmid_gff_output + "/{sample}/{gbk}/gff/recombinase_integrase_genes.gff",
    params:
        location = "{gbk}",
        seq_type = "plasmid"
    script:
        "scripts/python/DNAPlotter/find_interesting_gene_to_gff.py"



rule generate_template:
    input:
        gbk =plasmid_gbk + "/{sample}/{gbk}.gbk",
        orthologous_gene_track_template = rules.find_orthologous_gene.output[0],

    output:
        plasmid_gff_output + "/{sample}/{gbk}/DNAPlotter_template.txt"
    params:
        gff_dir = plasmid_gff_output + "/{sample}/{gbk}/gff",
        seq_type= "plasmid"
    conda: "coding_env"
    script:
        "scripts/python/DNAPlotter/generate_template.py"










