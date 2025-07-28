from io import StringIO
import os
from os import path

from Bio import AlignIO
import pandas as pd

def get_SNV_density():
    SNV_density_info_io = StringIO()
    SNV_density_info_io.write("GC_id\talignment_length_without_gap\tSNV_count\tSNV_density\n")

    for sample in fna_list:
        # Read the multiple sequence alignment
        alignment = AlignIO.read(sample, "fasta")

        GC_id = "_".join(path.basename(sample).split("_")[:2])

        # Array to hold the distinct columns
        distinct_columns = 0

        # Traverse through each column
        aln_len = alignment.get_alignment_length()
        for i in range(aln_len):
            column = alignment[:, i]
            distinct_elements = set(column)  # Convert the string to a set to get distinct elements
            # print(f"Column {i+1}: {distinct_elements} {len(distinct_elements)}")
            if len(distinct_elements) > 1:
                distinct_columns += 1

        SNV_density = distinct_columns / aln_len

        SNV_density_info_io.write(f"{GC_id}\t{aln_len}\t{distinct_columns}\t{SNV_density}\n")



    SNV_density_info_io.seek(0)
    df = pd.read_table(SNV_density_info_io, index_col=0)

    df.sort_values(by="SNV_density", ascending=False).to_csv(SNV_density_f, sep="\t")
    # 计算中位数（第二四分位数）
    median = df['SNV_density'].quantile(0.5)

    # 计算第一四分位数和第三四分位数
    q1 = df['SNV_density'].quantile(0.25)
    q3 = df['SNV_density'].quantile(0.75)
    print("SNV_density:")
    print("第一四分位数:", q1)
    print("中位数:", median)
    print("第三四分位数:", q3)

    print(df.loc[["GC_1826", "GC_1822", "GC_2202", "GC_1540", "GC_1827", "GC_1880", "GC_1821", "GC_2630", "GC_1828",
                  "GC_1829"], :].to_csv(sep="\t"))

    return df["SNV_density"]


def get_SAAV_density():
    SNV_density_info_io = StringIO()
    SNV_density_info_io.write("GC_id\talignment_length_without_gap\tSAAV_count\tSAAV_density\n")

    for sample in faa_list:
        # Read the multiple sequence alignment
        alignment = AlignIO.read(sample, "fasta")

        GC_id = "_".join(path.basename(sample).split("_")[:2])

        # Array to hold the distinct columns
        distinct_columns = 0
        aln_len_without_gap = 0

        # Traverse through each column
        aln_len = alignment.get_alignment_length()
        for i in range(aln_len):
            column = alignment[:, i]
            distinct_elements = set(column)  # Convert the string to a set to get distinct elements
            if "-" in distinct_elements:
                continue
            if len(distinct_elements) > 1:
                distinct_columns += 1
            aln_len_without_gap += 1

        SAAV_density = distinct_columns / aln_len_without_gap

        SNV_density_info_io.write(f"{GC_id}\t{aln_len_without_gap}\t{distinct_columns}\t{SAAV_density}\n")


    SNV_density_info_io.seek(0)
    df = pd.read_table(SNV_density_info_io, index_col=0)

    single_copy_gene_info_df = pd.read_table(single_copy_gene_info_file, index_col=1)

    df["bakta_Product"] = single_copy_gene_info_df.loc[df.index,"bakta_Product"]

    df.sort_values(by="SAAV_density", ascending=False).to_csv(SAAV_density_f, sep="\t")

    # 计算中位数（第二四分位数）
    median = df['SAAV_density'].quantile(0.5)

    # 计算第一四分位数和第三四分位数
    q1 = df['SAAV_density'].quantile(0.25)
    q3 = df['SAAV_density'].quantile(0.75)
    print("SAAV_density:")
    print("第一四分位数:", q1)
    print("中位数:", median)
    print("第三四分位数:", q3)

    print(df.loc[["GC_1826", "GC_1822", "GC_2202", "GC_1540", "GC_1827", "GC_1880", "GC_1821", "GC_2630", "GC_1828",
                  "GC_1829"], :].to_csv(sep="\t"))

    return df["SAAV_density"]




if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        sys.stdout = f
        sys.stderr = f


        fna_list = snakemake.input.fna_list
        faa_list = snakemake.input.faa_list
        single_copy_gene_info_file = snakemake.input.single_copy_gene_info_file


        SNV_density_f = snakemake.output.SNV_density
        SAAV_density_f = snakemake.output.SAAV_density
        merge_density_f = snakemake.output.merge_density

        pd.concat([get_SNV_density(),get_SAAV_density()],axis=1).to_csv(merge_density_f,sep="\t")


