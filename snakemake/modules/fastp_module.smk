"""
config key
    - reads_dir
    - unqualified
    - read_length
    - env_name
"""

FASTP_OUT_PATH = "fastp_output"
reads_dir = config["reads_dir"]
env_name = config["env_name"]
is_phred64 = config.get("is_phred64", False)
ffq_suffix = config.get("ffq_suffix","_1.fq.gz")
rfq_suffix = config.get("rfq_suffix","_2.fq.gz")
convert_phred33_output = "convert_phred33_output"


rule phred64_to_phred33:
    input:
        r1 = reads_dir + "/{sample}/{sample}" + ffq_suffix,
        r2 = reads_dir + "/{sample}/{sample}" + rfq_suffix,
    output:
        r1 = convert_phred33_output + "/{sample}/{sample}" + ffq_suffix,
        r2 = convert_phred33_output + "/{sample}/{sample}" + rfq_suffix,
    conda:
        env_name
    shell:
        "seqtk seq -Q 64 -V {input.r1} | gzip > {output.r1};"
        "seqtk seq -Q 64 -V {input.r2} | gzip > {output.r2};"



rule fastp:
    input:
        r1 = lambda wildcards: f"{convert_phred33_output if is_phred64 else reads_dir}/{wildcards.sample}/{wildcards.sample}{ffq_suffix}",
        r2 = lambda wildcards: f"{convert_phred33_output if is_phred64 else reads_dir}/{wildcards.sample}/{wildcards.sample}{rfq_suffix}",
    output:
        r1 = FASTP_OUT_PATH + "/{sample}/{sample}_paired_1.fq.gz",
        r2 = FASTP_OUT_PATH + "/{sample}/{sample}_paired_2.fq.gz",
        fo = FASTP_OUT_PATH + "/{sample}/{sample}_failed.fq.gz",
        json= ensure(temp(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.json"), non_empty=True),
        html = ensure(FASTP_OUT_PATH + "/{sample}/{sample}_fastp.html", non_empty=True), # Often, it is a good idea to combine ensure annotations with retry definitions, e.g. for retrying upon invalid checksums or empty files.
    log:
        "logs/fastp/{sample}_fastp.log"
    threads: 16
    params:
        unqualified=config['unqualified'],
        read_length=config['read_length'],
        N_base=config['N_base'] if 'N_base' in config else 3,
    conda:
        env_name
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --failed_out {output.fo} "
        "-j {output.json} -h {output.html} -w {threads} "
        "-q 20 " # the quality value that a base is qualified
        "-u {params.unqualified} " # how many percents of bases are allowed to be unqualified (0~100).
        "-n {params.N_base} " # number of N allowed
        "-l {params.read_length} " # reads shorter than length_required will be discarded, particularly those with adapter.
        "-5 " # enable trimming in 5' ends.
        "-3 " # enable trimming in 3' ends.
        "&> {log}"


