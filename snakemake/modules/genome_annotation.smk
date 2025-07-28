"""
config key
    - genome_dir
    - env_name
    - samples
    - busco_database
    - checkM_output
    - busco_output
"""


genome_suffix = config.get("genome_suffix","fasta")

bakta_output = config["bakta_output"]
pgap_output = config["pgap_output"]

bakta_db = config["bakta_db"]
genome_dir = config["genome_dir"]
bakta_replicons_table_dir = config["bakta_replicons_table_dir"]
bakta_proteins_output = config["bakta_proteins_output"]
env_name = config["env_name"]
SAMPLES = config["samples"]



# Check whether the environment variable of the pgap database path is configured
envvars:
    "PGAP_INPUT_DIR"



rule pgap_yaml:
    input:
        genome=genome_dir + "/{sample}." + genome_suffix,
    output:
        generic=pgap_output + "/{sample}/generic.yaml",
        metadata=pgap_output + "/{sample}/metadata.yaml",
        genome=temp(pgap_output + "/{sample}/{sample}.fasta"),
    params:
        genus_species=lambda wildcards: config["organism"][wildcards.sample],
    run:
        
        shell(f"cp {input.genome} {output.genome}")
        generic_yaml = f"""fasta:
    class: File
    location: {wildcards.sample}.fasta
submol:
    class: File
    location: metadata.yaml"""

        metadata_yaml = f"""organism:
    genus_species: {params["genus_species"]}
    strain: {wildcards.sample}
locus_tag_prefix: {wildcards.sample}
authors:
    - author:
        first_name: 'Jianmin'
        last_name: 'Xie'
contact_info:
    first_name: 'Jianmin'
    last_name: 'Xie'
    email: '19jmxie@stu.edu.cn'
    organization: 'Shantou University'
    department: 'College of Science'
    street: '243 Daxue Road'
    city: 'Shantou'
    state: 'Guangdong'
    postal_code: '515063'
    country: 'China'"""
        with open(output["generic"],"w") as f:
            f.write(generic_yaml)
        with open(output["metadata"],"w") as f:
            f.write(metadata_yaml)





rule pgap:
    input:
        generic = rules.pgap_yaml.output.generic,
        genome = rules.pgap_yaml.output.genome,
    output:
        gff = pgap_output + "/{sample}/result/annot.gff",
        gbk = pgap_output + "/{sample}/result/annot.gbk",
        faa = pgap_output + "/{sample}/result/annot.faa",
        translated_faa = pgap_output + "/{sample}/result/annot_translated_cds.faa"
    params:
        threads=128,
        out_dir = pgap_output + "/{sample}/result"
    benchmark:
        "benchmarks/pgap/{sample}_benchmark.tsv"
    threads: workflow.cores / 4
    log:
        "logs/pgap/{sample}_pgap.log"
    shell:
        "rm -rf {params.out_dir};"
        "pgap.py --no-self-update --taxcheck --ignore-all-errors -c {params.threads} -n -o {params.out_dir} {input.generic} &> {log}"





rule bakta_proteins:
    input:
        rules.pgap.output.faa
    output:
        bakta_proteins_output + "/{sample}/{sample}.tsv"
    params:
        out_dir=bakta_proteins_output + "/{sample}",
        db=bakta_db,
    log:
        "logs/bakta_proteins/{sample}_bakta_proteins.log"
    threads: workflow.cores / 4
    conda:
        env_name
    shell:
        "bakta_proteins -f "
        "--db {params.db} "
        "--verbose "
        "--output {params.out_dir} "
        "--prefix {wildcards.sample} "
        "--threads {threads} {input} "
        "&> {log}"








