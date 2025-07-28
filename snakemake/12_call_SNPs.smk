"""
name: call SNPs with reference genome
date: 2024-04-07
description: Snakemake pipeline to call SNPs by Taking all genomes in turn as reference genomes
author: Jianmin Xie
dependencies:
    - bowtie2 = 2.5.3
    - freebayes = 1.3.6 # https://github.com/freebayes/freebayes
    - snpEff = 5.2
    - pandas = 2.2.3
rules:
    -
"""
import os
from os import path
import pandas as pd
shell.executable('/bin/zsh') # set a shell to run commands

configfile: "configs/call_SNPs_config.yaml"

env_name = config["env_name"]
genome_dir = config["genome_dir"]
gbk_dir = config["gbk_dir"]
single_copy_gene_ann_dir = config["single_copy_gene_ann_dir"]
fastp_output = config["fastp_output"]
is_phred64 = config.get("is_phred64",False)

bowtie2_index = "bowtie2_index"
bowtie2_output = "bowtie2_output"
call_variants_output = "call_variants_output"
variants_annotation_output = "variants_annotation_output"
variants_annotation_summary_output = "variants_annotation_summary"

blast_db_output = "blast_db"



SAMPLES =  sorted([".".join(i.split(".")[:-1]) for i in os.listdir(genome_dir) if i.endswith("fasta")])
# print(SAMPLES)


def get_target_res(w):
    t_list = []

    for sample in SAMPLES:
        t_list += expand(bowtie2_output + f"/{sample}" + "/{ref_sample}/{ref_sample}_idxstats.tsv", ref_sample=[i for i in SAMPLES if i != sample])

    return t_list

rule all:
    input:
        # expand(bowtie2_index + "/{sample}/{sample}.1.bt2", sample=SAMPLES),
        get_target_res,
        expand(variants_annotation_output + "/config/data/{sample}/snpEffectPredictor.bin", sample=SAMPLES),
        expand(variants_annotation_summary_output + "/all/{type}.tsv", type=["average_SNV_count", "average_SNV_density_Kbp", "AA_changes_percentage"]),
    default_target: True


# Ref：https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
rule build_bowtie2_index:
    input:
        genome_dir + "/{sample}.fasta"
    output:
        bowtie2_index + "/{sample}/{sample}.1.bt2"
    params:
        index_prifix=bowtie2_index + "/{sample}/{sample}"
    threads: 32
    log:
        "logs/build_bowtie2_index/{sample}.log"
    conda:
        env_name
    shell:
        "bowtie2-build --threads {threads} {input} {params.index_prifix} &> {log}"

# Ref：https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
rule bowtie2:
    input:
        r1=fastp_output + "/{sample}/{sample}_paired_1.fq.gz",
        r2=fastp_output + "/{sample}/{sample}_paired_2.fq.gz",
        bowtie2_index = bowtie2_index + "/{ref_sample}/{ref_sample}.1.bt2"
    output:
        temp(bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}.sam")
    params:
        index_prifix=bowtie2_index + "/{ref_sample}/{ref_sample}",
        is_phred64 = "--phred64" if is_phred64 else ""
    log:
        "logs/bowtie2/{sample}/{ref_sample}_bowtie2.log"
    conda:
        env_name
    threads: 64
    shell:
        "bowtie2 -x {params.index_prifix} -1 {input.r1} -2 {input.r2} -S {output} -p {threads} {params.is_phred64} &> {log}"

# Sorting and compression
rule samtools_sort:
    input:
        rules.bowtie2.output
    output:
        bam = temp(bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}.sort.bam"),
        bai = temp(bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}.sort.bam.bai"),

    threads: 24
    conda:
        env_name
    log:
        "logs/samtools_sort/{sample}/{ref_sample}_samtools.log"
    shell:
        "samtools sort -O BAM -o {output.bam} -@ {threads} {input} &> {log};"
        "samtools index -@ {threads} {output.bam};"

rule mapping_summary:
    input:
        bam = rules.samtools_sort.output.bam,
        bai = rules.samtools_sort.output.bai,
    output:
        idxstats = bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}_idxstats.tsv",
        flagstat = bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}_flagstat.txt",
        coverage = bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}_coverage.tsv",
        coverage_histogram = bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}_coverage_histogram.tsv",
        depth_histogram = bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}_depth_histogram.tsv",
    conda:
        env_name
    shell:
        "samtools idxstats {input.bam} > {output.idxstats};"
        "samtools flagstat {input.bam} > {output.flagstat};"
        "samtools coverage {input.bam} > {output.coverage};"
        "samtools coverage -A {input.bam} > {output.coverage_histogram};"
        "samtools coverage -A --plot-depth {input.bam} > {output.depth_histogram};"


# Filter out the unmatched

rule samtools_filter:
    input:
        rules.samtools_sort.output.bam
    output:
        bowtie2_output + "/{sample}/{ref_sample}/{ref_sample}.mapped.sort.bam"

    threads: 24
    conda:
        env_name
    log:
        "logs/samtools_filter/{sample}/{ref_sample}_samtools.log"
    shell:
        "samtools view -F 4 -@ {threads} -bS -o {output} {input} &> {log};" # -F 4 滤掉没有mapping上records
        "samtools index -@ {threads} {output};"



rule freebayes:
    input:
        bam = rules.samtools_filter.output,
        ref_genome = genome_dir + "/{ref_sample}.fasta"
    output:
        call_variants_output + "/{sample}/{ref_sample}.vcf"
    params:
        min_coverage = 10, # Mapping Depth
        min_base_quality = 20
        min_alternate_fraction = 0.75, # The minimum ratio of the mutated bases in this series to all bases
        haplotype_length = -1 # Do not test haplotypes for consecutive mutations
    conda:
        env_name
    shell:
        "freebayes -f {input.ref_genome} --min-coverage {params.min_coverage} --min-base-quality {params.min_base_quality} --min-alternate-fraction {params.min_alternate_fraction} --haplotype-length {params.haplotype_length} --vcf {output}  {input.bam}"



rule snpEff_config:
    input:
        ref_gbk = gbk_dir + "/{sample}/result/annot.gbk",
        ref_genome = genome_dir + "/{sample}.fasta"
    output:
        config_file = variants_annotation_output + "/config/{sample}_snpEff.config",
        ref_gbk = variants_annotation_output + "/config/data/{sample}/genes.gbk"
    run:
        shell("cp {input.ref_gbk} {output.ref_gbk}")
        base_str = f"""codon.Bacterial_and_Plant_Plastid : TTT/F , TTC/F , TTA/L , TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G

{wildcards.sample}.genome : {wildcards.sample}
"""
        for line in shell("grep '>' {input.ref_genome}",iterable=True):
            base_str += f"{wildcards.sample}.{line.strip().split(' ')[0][1:]}.codonTable : Bacterial_and_Plant_Plastid\n"
        with open(output.config_file, "w") as f:
            f.write(base_str)


rule build_snpEff_database:
    input:
        rules.snpEff_config.output.config_file
    output:
        variants_annotation_output + "/config/data/{sample}/snpEffectPredictor.bin"
    conda:
        env_name
    log:
        "logs/build_snpEff_database/{sample}.log"
    shell:
        "snpEff build -genbank -c {input} -v {wildcards.sample} &> {log}"

rule run_snpEff:
    input:
        config_file = variants_annotation_output + "/config/{ref_sample}_snpEff.config",
        vcf_file = rules.freebayes.output
    output:
        ann_vcf = variants_annotation_output + "/result/{sample}/{ref_sample}/{ref_sample}_ann.vcf",
        summary_html = variants_annotation_output + "/result/{sample}/{ref_sample}/{ref_sample}_ann_summary.html"
    log:
        "logs/run_snpEff/{sample}/{ref_sample}_snpEff.log"
    conda:
        env_name
    shell:
        "(snpEff -no-downstream -no-upstream -no-utr -no-intergenic -c {input.config_file} -htmlStats {output.summary_html} -v {wildcards.ref_sample} {input.vcf_file} > {output.ann_vcf}) &> {log}"


rule summary_variants_annotation:
    input:
        ann_vcf_list = lambda w: expand(variants_annotation_output + f"/result/{w.sample}/{{ref_sample}}/{{ref_sample}}_ann.vcf", ref_sample=[i for i in SAMPLES if i != w.sample]),
    output:
        SNV_count_file = variants_annotation_summary_output + "/{sample}/SNV_count.tsv",
        SNV_density_file = variants_annotation_summary_output + "/{sample}/SNV_density.tsv",
        AA_changes_percentage_file = variants_annotation_summary_output + "/{sample}/AA_changes_percentage.tsv",
        # base_changes_info_dir = directory(variants_annotation_summary_output + "/base_changes_info"),
        # SNV_codon_position_info_dir = directory(variants_annotation_summary_output + "/SNV_codon_position_info"),
        # SNV_variant_type_info_dir = directory(variants_annotation_summary_output + "/SNV_variant_type_info"),
        # AA_changes_info_dir = directory(variants_annotation_summary_output + "/AA_changes_info"),
    params:
        single_copy_gene_ann_dir = single_copy_gene_ann_dir
    conda:
        env_name
    log:
        "logs/summary_variants_annotation/{sample}.log"
    script:
        "scripts/python/summary_variants_annotation_all_ref.py"

rule summary_all:
    input:
        SNV_count_file_list = expand(variants_annotation_summary_output + "/{sample}/SNV_count.tsv", sample=SAMPLES),
        SNV_density_file_list = expand(variants_annotation_summary_output + "/{sample}/SNV_density.tsv", sample=SAMPLES),
        AA_changes_percentage_file_list = expand(variants_annotation_summary_output + "/{sample}/AA_changes_percentage.tsv", sample=SAMPLES),

    output:
        average_SNV_count_file = variants_annotation_summary_output + "/all/average_SNV_count.tsv",
        average_SNV_density_file = variants_annotation_summary_output + "/all/average_SNV_density_Kbp.tsv",
        AA_changes_percentage_file = variants_annotation_summary_output + "/all/AA_changes_percentage.tsv",
        AA_changes_percentage_top25_file = variants_annotation_summary_output + "/all/AA_changes_percentage_top25.tsv",
    run:
        data_dict = {
            "SNV_count_file_list":{"column":"Average count", "out":"average_SNV_count_file"},
            "SNV_density_file_list":{"column":"Average density", "out":"average_SNV_density_file"},
            "AA_changes_percentage_file_list":{"column":"Percentage (%)", "out":"AA_changes_percentage_file"},
        }
        for key in data_dict:
            S_list = []
            for SNV_count_file in input[key]:
                S = pd.read_table(SNV_count_file,index_col=0)[data_dict[key]["column"]]
                S.name = SNV_count_file.split("/")[-2]
                S_list.append(S)
            S_df = pd.concat(S_list,axis=1).fillna(0)
            S_df = S_df.loc[S_df.mean(axis=1).sort_values(ascending=False).index]
            S_df.to_csv(output[data_dict[key]["out"]],sep="\t")
            if key == "AA_changes_percentage_file_list":
                S_df.T.iloc[:,:25].to_csv(output.AA_changes_percentage_top25_file,sep="\t")

