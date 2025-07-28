"""
config key
    - genome_dir
    - env_name
    - samples
    - busco_database
    - checkM_output
    - busco_output
"""

import os

# 默认值
genome_suffix = config.get("genome_suffix","fasta")

bakta_output = "bakta_output"

bakta_db = config["bakta_db"]
genome_dir = config["genome_dir"]
env_name = config["env_name"]
SAMPLES = [".".join(i.split(".")[:-1]) for i in os.listdir(genome_dir) if not i.startswith(".")]



# 检查有没有配置pgap数据库路径的环境变量
envvars:
    "PGAP_INPUT_DIR"

rule bakta:
    input:
        genome=genome_dir + "/{sample}." + genome_suffix
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
        "--threads {threads} {input.genome} "
        "&> {log}"




