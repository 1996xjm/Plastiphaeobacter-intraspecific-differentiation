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
shell.executable('/bin/zsh')

chromosome_gbk = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/genome/chromosome_gbk'
xmfa_file = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/Mauve_output/aln_res.xmfa'
backbone_file = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/Mauve_output/aln_res.backbone'

ppanggolin_GIs_file = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/ppanggolin/chromosome/regions_of_genomic_plasticity.tsv'

mcl_cluster_file = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/my_pangenome/pan_genome_output/table/mcl-clusters.tsv'
ref_annotation_dir = f'/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/merge_annotation'

chromosome_phage_regions_dir = '/Users/xjm/Documents/sampling/final_representative_strain/Phaeobacter_gallaeciensis/TG_genome/online_annotation/PHASTEST/chromosome'
SAMPLES = sorted([i.replace("_chromosome_1.gbk","") for i in os.listdir(chromosome_gbk) if i.endswith("gbk")])

chromosome_gff_output = "chromosome"
print(SAMPLES)


rule all:
    input:
        expand(chromosome_gff_output + "/{sample}/gff/GIs.gff", sample=SAMPLES),
        expand(chromosome_gff_output + "/{sample}/orthologous_gene_track_template.tsv", sample=SAMPLES),
        expand(chromosome_gff_output + "/{sample}/gff/{type}.gff", sample=SAMPLES, type=["transposase_genes","recombinase_integrase_genes", "chromosome_phage_regions_intact","chromosome_phage_regions_questionable","chromosome_phage_regions_incomplete"]),
        expand(chromosome_gff_output + "/{sample}/DNAPlotter_template.txt", sample=SAMPLES),


rule extract_genomic_islands:
    input:
        ppanggolin_GIs_file = ppanggolin_GIs_file,
    output:
        chromosome_gff_output + "/{sample}/gff/GIs.gff"
    script:
        "scripts/python/DNAPlotter/extract_ppanggolin_GIs_to_gff.py"



rule find_orthologous_gene:
    input:
        mcl_cluster_file = mcl_cluster_file,
        ref_annotation_file = ref_annotation_dir + "/{sample}_annotation.tsv",


    output:
        chromosome_gff_output + "/{sample}/orthologous_gene_track_template.tsv",
    params:
        all_smaples = SAMPLES,
        gff_out_dir = chromosome_gff_output + "/{sample}/gff",
        track_start_position = 0.7,
        location= "{sample}_chromosome_1"
    script:
        "scripts/python/DNAPlotter/find_orthologous_gene_to_gff.py"



rule find_interesting_gene:
    input:
        ref_annotation_file = ref_annotation_dir + "/{sample}_annotation.tsv",
        chromosome_phage_regions = chromosome_phage_regions_dir + "/{sample}/predicted_phage_regions.json",

    output:
        transposase_gene_gff_file = chromosome_gff_output + "/{sample}/gff/transposase_genes.gff",
        recombinase_integrase_gene_gff_file = chromosome_gff_output + "/{sample}/gff/recombinase_integrase_genes.gff",
        chromosome_GTA_gff_file = chromosome_gff_output + "/{sample}/gff/chromosome_GTA.gff",
        chromosome_phage_regions_intact_gff_file = chromosome_gff_output + "/{sample}/gff/chromosome_phage_regions_intact.gff",
        chromosome_phage_regions_questionable_gff_file = chromosome_gff_output + "/{sample}/gff/chromosome_phage_regions_questionable.gff",
        chromosome_phage_regions_incomplete_gff_file = chromosome_gff_output + "/{sample}/gff/chromosome_phage_regions_incomplete.gff",
    params:
        location = "{sample}_chromosome_1",
        seq_type= "chromosome"
    script:
        "scripts/python/DNAPlotter/find_interesting_gene_to_gff.py"



rule generate_template:
    input:
        gbk =chromosome_gbk + "/{sample}_chromosome_1.gbk",
        orthologous_gene_track_template = rules.find_orthologous_gene.output[0],

    output:
        chromosome_gff_output + "/{sample}/DNAPlotter_template.txt"
    params:
        gff_dir = chromosome_gff_output + "/{sample}/gff",
        seq_type= "chromosome"
    conda: "coding_env"
    script:
        "scripts/python/DNAPlotter/generate_template.py"










