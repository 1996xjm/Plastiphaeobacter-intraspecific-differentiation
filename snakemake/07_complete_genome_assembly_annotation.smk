"""
name: Third generation sequence assembly
date: 2023-04-24
description: Snakemake pipeline to do Third generation sequence assembly
author: Jianmin Xie
dependencies:
    - fastp = 0.23.2
    - flye = 2.9.2
    - minimap2 = 2.24
    - samtools = 1.17
    - pilon = 1.24
    - seqtk = 1.4 # convert phred64 to 33
    - pyyaml = 6.0.1
    - quast = 5.2.0
    - pgap = v2024-04-27.build7426
    - bakta = 1.9.2
    - prodigal = 2.6.3
    - biopython = 1.83
    - blast = 2.15.0
    - pandas = 2.2.1

"""

import os
from os import path
configfile: "configs/TG_assembly_config.yaml"

pilon_jar = config["pilon_jar"]
SG_reads_dir = config["SG_reads_dir"]
TG_reads_dir = config["TG_reads_dir"]
polish_times = config["polish_times"] # polishing recurse time
env_name = config["env_name"]
bakta_db = config["bakta_db"]
SAMPLES = os.listdir(TG_reads_dir)

FASTP_OUT_PATH = "fastp_output"
flye_OUT_PATH = "flye_output"
minimap2_OUT_PATH = "minimap2_output"
pilon_OUT_PATH = "pilon_output"
final_genome_PATH = "genome/complete"
final_plasmid_PATH = "genome/plasmid"
final_chromosome_PATH = "genome/chromosome"
final_chromosome_gbk_PATH = "genome/chromosome_gbk"
final_plasmid_gbk_PATH = "genome/plasmid_gbk"
final_plasmid_faa_PATH = "genome/plasmid_faa"
unfixed_genome_PATH = "unfixed_genome"
bakta_replicons_table_dir = "bakta_replicons_table"
checkM_output = "quality_assessment/checkM_output"
busco_output = "quality_assessment/busco_output"
quast_genome_output = "quality_assessment/quast_genome_output"
quast_chromosome_output = "quality_assessment/quast_chromosome_output"
bakta_output = "bakta_for_adjust_genome"
pgap_output = "pgap_output"
bakta_proteins_output = "bakta_proteins_output"


shell.executable('/bin/zsh') # set a shell to run commands

module fastp_module:
    # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
    snakefile: "modules/fastp_module.smk"
    config: config["fastp_config"]

quality_assessment_config = {
    "genome_dir": final_genome_PATH,
    "chromosome_dir": final_chromosome_PATH,
    "env_name": config["quality_assessment_env_name"],
    "samples": SAMPLES,
    "busco_database": config["busco_database"],
    "checkM_output": checkM_output,
    "busco_output": busco_output,
    "quast_genome_output": quast_genome_output,
    "quast_chromosome_output": quast_chromosome_output,

}

module quality_assessment_m:
    snakefile: "modules/quality_assessment_module.smk"
    config: quality_assessment_config

genome_annotation_config = {
        "bakta_db":bakta_db,
        "genome_dir":final_genome_PATH,
        "bakta_replicons_table_dir":bakta_replicons_table_dir,
        "organism":config["organism"],
        "env_name":env_name,
        "samples": SAMPLES,
        "bakta_output": bakta_output,
        "pgap_output": pgap_output,
        "bakta_proteins_output": bakta_proteins_output,
    }
module genome_annotation_m:
    snakefile: "modules/genome_annotation.smk"
    config: genome_annotation_config

# use fastp module
use rule * from fastp_module as my_*
use rule * from quality_assessment_m as qa_*
use rule * from genome_annotation_m as ga_*

rule all:
    input:
        busco_output + "/batch_summary.tsv",
        checkM_output + "/checkm_result.tsv",
        quast_genome_output + "/report.tsv",
        quast_chromosome_output + "/report.tsv",
        # expand(bakta_output + "/{sample}/{sample}.tsv",sample=SAMPLES),
        # expand(pgap_output + "/{sample}/result/annot.tsv",sample=SAMPLES),
        # expand(final_genome_PATH + "/{sample}.fasta",sample=SAMPLES),
        expand(bakta_proteins_output + "/{sample}/{sample}.tsv",sample=SAMPLES),
        expand(pgap_output + "/{sample}/result/annot.gff",sample=SAMPLES),
        final_plasmid_PATH,
        final_plasmid_faa_PATH,
        final_chromosome_PATH,
        final_chromosome_gbk_PATH,
        final_plasmid_gbk_PATH
    default_target: True








# Mapping
rule flye:
    input:
        TG_reads_dir + "/{sample}/{sample}.pass.fastq.gz"
    output:
        ensure(flye_OUT_PATH + "/{sample}/assembly.fasta", non_empty=True)
    params:
        genomeSize = lambda wildcards: config["genome_size"][wildcards.sample],
        outdir = flye_OUT_PATH + "/{sample}",
        asm_coverage = 300 # Using longest 300x reads for contig assembly
    threads: workflow.cores * (0.3 if len(SAMPLES) > 1 else 1)
    log:
        "logs/flye/{sample}_flye.log"
    benchmark:
        "benchmarks/flye/{sample}_benchmark.tsv"
    conda:
        env_name
    shell:
        "flye --nano-raw {input} -g {params.genomeSize} --asm-coverage {params.asm_coverage} --out-dir {params.outdir} --threads {threads} &> {log}"




# Correct the assembly


def recurse_genome(wcs):
    n = int(wcs.n)
    if n == 1:
        return flye_OUT_PATH + f"/{wcs.sample}/assembly.fasta"
    elif n > 1:
        return pilon_OUT_PATH + f"/time_{n-1}/{wcs.sample}/{wcs.sample}_pilon_{n-1}.fasta"
    else:
        raise ValueError(f"loop numbers must be 1 or greater: received {wcs.n}")



# Polishing

rule minimap2:
    input:
        genome = recurse_genome,
        r1 = rules.my_fastp.output.r1,
        r2 = rules.my_fastp.output.r2,
    output:
        minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_aln.bam"
    threads: 10
    log:
        "logs/minimap2/time_{n}/{sample}_minimap2.log"
    conda:
        env_name
    shell:
        "(minimap2 -ax sr {input.genome} {input.r1} {input.r2} | samtools sort -@ {threads} -O bam -o {output}) &> {log};"
        "samtools index -@ {threads} {output};"

rule mapping_summary:
    input:
        rules.minimap2.output
    output:
        idxstats = minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_idxstats.tsv",
        flagstat = minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_flagstat.txt",
        coverage = minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_coverage.tsv",
        coverage_histogram = minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_coverage_histogram.txt",
        depth_histogram = minimap2_OUT_PATH + "/time_{n}/{sample}/{sample}_depth_histogram.txt",
    conda:
        env_name
    shell:
        "samtools idxstats {input} > {output.idxstats};"
        "samtools flagstat {input} > {output.flagstat};"
        "samtools coverage {input} > {output.coverage};"
        "samtools coverage -A {input} > {output.coverage_histogram};"
        "samtools coverage -A --plot-depth {input} > {output.depth_histogram};"

rule pilon:
    input:
        genome = recurse_genome,
        bam = rules.minimap2.output,
        summary = rules.mapping_summary.output,
    output:
        pilon_OUT_PATH + "/time_{n}/{sample}/{sample}_pilon_{n}.fasta"
    params:
        prefix = "{sample}_pilon_{n}",
        outdir= pilon_OUT_PATH + "/time_{n}/{sample}",
        pilonJar = pilon_jar
    log:
        "logs/pilon/time_{n}/{sample}_pilon.log"
    shell:
        "java -Xmx16G -jar {params.pilonJar} --genome {input.genome} --frags {input.bam} --output {params.prefix} --fix all --mindepth 10 --minqual 20 --outdir {params.outdir} --changes --verbose &> {log}"



"""
# Parameters Description :
--mindepth depth
          Variants (snps and indels) will only be called if there is coverage of good pairs
          at this depth or more; if this value is >= 1, it is an absolute depth, if it is a
          fraction < 1, then minimum depth is computed by multiplying this value by the mean
          coverage for the region, with a minumum value of 5 (default 0.1: min depth to call
          is 10% of mean coverage or 5, whichever is greater).
--minqual
              Minimum base quality to consider for pileups (default 0)
"""







# Rename contigs
rule rename:
    input:
        pilon_OUT_PATH + f"/time_{polish_times}/{{sample}}/{{sample}}_pilon_{polish_times}.fasta"
    output:
        fasta=unfixed_genome_PATH + "/{sample}_unfixed.fasta",
        bakta_replicons_table=bakta_replicons_table_dir + "/{sample}_replicons.tsv"
    params:
        sample = "{sample}",
        assembly_info = flye_OUT_PATH + "/{sample}/assembly_info.txt",
        polish_times = polish_times
    log:
        "logs/rename/{sample}_rename.log"
    conda:
        env_name
    script:
        "scripts/python/rename_TG_genome.py"






# Identify the locus of gene dnaA
rule bakta:
    input:
        genome=rules.rename.output.fasta,
        bakta_replicons_table=bakta_replicons_table_dir + "/{sample}_replicons.tsv"
    output:
        bakta_output + "/{sample}/{sample}.tsv"
    params:
        out_dir=bakta_output + "/{sample}",
        db=bakta_db,
        genus=lambda wildcards: config["organism"][wildcards.sample].split(" ")[0],
        species=lambda wildcards: config["organism"][wildcards.sample].split(" ")[1],
    log:
        "logs/bakta/{sample}_bakta.log"
    threads: workflow.cores / 4
    conda:
        env_name
    shell:
        "bakta -f "
        "--db {params.db} "
        "--verbose "
        "--genus {params.genus} "
        "--species {params.species} "
        "--output {params.out_dir} "
        "--prefix {wildcards.sample} "
        "--locus-tag {wildcards.sample} "
        "--replicons {input.bakta_replicons_table} "
        "--threads {threads} {input.genome} "
        "&> {log}"


# Adjust the fasta file

rule adjust_genome:
    input:
        genome = rules.rename.output.fasta,
        gff = rules.bakta.output
    output:
        final_genome_PATH + "/{sample}.fasta"
    log:
        "logs/adjust_genome/{sample}_adjust_genome.log"
    conda:
        env_name
    script:
        "scripts/python/adjust_TG_genome_by_dnaA.py"

# Extract plasmid and chromosome to separate fasta files

rule extract_plasmid_chromosome:
    input:
        expand(final_genome_PATH + "/{sample}.fasta", sample=SAMPLES)
    output:
        plasmid_dir = directory(final_plasmid_PATH),
        chromosome_dir = directory(final_chromosome_PATH),
        chromosome_fasta = expand(final_chromosome_PATH + "/{sample}_chromosome_1.fasta", sample=SAMPLES),
    conda:
        env_name

    script:
        "scripts/python/extract_plasmid_chromosome.py"


rule extract_plasmid_faa:
    input:
        translated_faa = expand(pgap_output + "/{sample}/result/annot_translated_cds.faa",sample=SAMPLES),
        annot_gbk = expand(pgap_output + "/{sample}/result/annot.gbk",sample=SAMPLES),
    output:
        faa = directory(final_plasmid_faa_PATH),
        chromosome_gbk = directory(final_chromosome_gbk_PATH),
        plasmid_gbk = directory(final_plasmid_gbk_PATH),
    conda:
        env_name

    script:
        "scripts/python/extract_plasmid_faa.py"






