from os import path
import os
from Bio import SeqIO
import datetime
from log import Log

def filterLength(fastaPath, outputPath, logPath, thLength=1000, thCoverage=10):
    """
    滤掉长度低于某个阈值的contigs
    1000是细菌单个基因的平均长度
    覆盖度低于10，说明这条contig很有可能是污染序列，或者这个地方确实是测序深度不够
    上传基因组到NCBI会帮你检测污染序列，PT34_32被鉴定到2条contigs是污染序列，其覆盖度低于10，其他正常contigs都高于10
    """
    log_f = open(logPath, "w")
    genomeDir = path.dirname(outputPath)
    if not path.exists(genomeDir):
        os.mkdir(genomeDir)

    if not path.exists(fastaPath):
        print("No such fasta file")
        print(fastaPath)
        os._exit(1)
    sampleName = ".".join(path.basename(outputPath).split(".")[:-1])
    records = SeqIO.parse(fastaPath, "fasta")
    good_records = []

    total_count = 0
    total_length = 0

    success_count = 0
    success_length = 0
    for rec in records:
        seqLength = len(rec.seq)
        seqIDSplit = rec.id.split('_')
        seqCoverage = float(seqIDSplit[-1])
        total_count += 1
        total_length += seqLength

        if seqLength >= thLength and seqCoverage >= thCoverage:
            rec.description = rec.id
            rec.id = f"{sampleName}_contig_{seqIDSplit[1]}"
            good_records.append(rec)
            success_length += seqLength
            success_count += 1
        else:
            Log.warning(f"{rec.id} fail to pass.", log_f)
    with open(outputPath, "w") as f:
        SeqIO.write(good_records, f, "fasta")

    out_res = f"""
    ## Before filter:
    total number of contigs: {total_count}
    total length of contigs: {total_length} bp

    ## After filter:
    total number of contigs: {success_count}
    total length of contigs: {success_length} bp
    """
    print(out_res, file=log_f)
    log_f.close()


if __name__ == '__main__':
    filterLength(snakemake.input[0], snakemake.output[0], snakemake.log[0], snakemake.params.thLength,
                 snakemake.params.thCoverage)