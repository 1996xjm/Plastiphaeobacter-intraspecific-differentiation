"""
config key
    - genome_dir
    - env_name
    - samples
    - busco_database
    - checkM_output
    - busco_output
"""

genome_dir = config["genome_dir"]
chromosome_dir = config["chromosome_dir"]
genome_suffix = config.get("genome_suffix","fasta")
checkM_output_path = config["checkM_output"]
quast_genome_output_path = config["quast_genome_output"]
quast_chromosome_output_path = config["quast_chromosome_output"]
busco_output_path = config["busco_output"]
busco_database = config.get("busco_database", "bacteria_odb10")
env_name = config["env_name"]
SAMPLES = config["samples"]




rule checkM:
    input:
        expand(genome_dir + "/{sample}." + genome_suffix, sample=SAMPLES)
    output:
        res=ensure(checkM_output_path + "/checkm_result.tsv", non_empty=True),
        tem=temp(directory(checkM_output_path + "/")),
    params:
        dir = genome_dir,
        out = checkM_output_path,
        suffix = genome_suffix
    threads: len(SAMPLES) if len(SAMPLES) < workflow.cores else workflow.cores
    conda:
        env_name
    log:
        "logs/checkM/checkM.log"
    shell:
        "checkm lineage_wf -f {output.res} --tab_table -x {params.suffix} -t {threads} --pplacer_threads {threads} {params.dir} {output.tem} &> {log}"


rule busco:
    input:
        expand(genome_dir + "/{sample}." + genome_suffix, sample=SAMPLES)
    output:
        ensure(busco_output_path + "/batch_summary.txt", non_empty=True)
    params:
        dir = genome_dir,
        out = busco_output_path,
        database = busco_database,
        download_path = busco_output_path + "/busco_downloads"
    conda:
        env_name
    log:
        "logs/busco/busco.log"
    shell:
        "busco -i {params.dir}  -l {params.database} --download_path {params.download_path} -m genome -o {params.out} -f &> {log}"


rule rename_busco:
    input:
        rules.busco.output
    output:
        ensure(busco_output_path + "/batch_summary.tsv",non_empty=True)
    shell:
        "mv {input} {output}"



rule quast_genome:
    input:
        expand(genome_dir + "/{sample}." + genome_suffix, sample=SAMPLES)
    output:
        ensure(quast_genome_output_path + "/report.tsv", non_empty=True)
    params:
        out = quast_genome_output_path
    conda:
        "quality_assessment"
    shell:
        "quast.py {input} -o {params.out} &> /dev/null"



rule quast_chromosome:
    input:
        expand(chromosome_dir + "/{sample}_chromosome_1." + genome_suffix, sample=SAMPLES)
    output:
        ensure(quast_chromosome_output_path + "/report.tsv", non_empty=True)
    params:
        out = quast_chromosome_output_path
    conda:
        "quality_assessment"
    shell:
        "quast.py {input} -o {params.out} &> /dev/null"




