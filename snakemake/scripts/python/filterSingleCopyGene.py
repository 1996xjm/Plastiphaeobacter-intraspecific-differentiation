from os import path
import os
from Bio import SeqIO
import datetime
from log import Log

def filterLength(single_copy_gene_dir, outputPath, logPath):
    log_f = open(logPath, "w")
    Log.info(f"{'#' * 10}Filter Single Copy Gene{'#' * 10}\n{'-' * 60}", log_f)
    OG_list = [i for i in os.listdir(single_copy_gene_dir) if not i.startswith(".")]
    Log.info(f"Total {len(OG_list)} single copy OGs.", log_f)
    count = 0  # 用于存放序列有差异的OG



    # 把有差异的OG复制到另一个文件夹
    if not path.exists(outputPath):
        os.mkdir(outputPath)

    for OG in OG_list:
        # print(OG)
        record_list = list(SeqIO.parse(path.join(single_copy_gene_dir, OG), "fasta"))  # 迭代对象要转换成list才能重复使用
        seq_record = [i.seq for i in record_list]
        flag = False
        # 判断所有序列是否有差异，有才用于串联建树
        for i in seq_record:
            if seq_record[0] != i:
                flag = True
                count += 1
                break

        if flag:
            """去序列名称后面的数字"""
            for rec in record_list:
                new_id = "_".join(rec.id.split("_")[:-1])
                rec.id = new_id
                rec.description = ""
            SeqIO.write(record_list, path.join(outputPath, OG), "fasta")

        if not flag:
            Log.info(f"Identital single copy OGs: {OG}", log_f)

    Log.info(f"Remove {len(OG_list) - count} identital single copy OGs.", log_f)
    log_f.close()


if __name__ == '__main__':
    filterLength(snakemake.params.single_copy_gene_dir, snakemake.output[0], snakemake.log[0])