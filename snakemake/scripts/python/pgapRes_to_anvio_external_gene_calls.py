
import pandas as pd
import numpy as np
from Bio import SeqIO


def faa_to_dict(faa_file):
    record_dict = SeqIO.index(faa_file, "fasta")

    faa_dict = {}

    for key in record_dict.keys():
        new_key = "_".join(key.split("_prot_")[1].split("_")[:-1])
        faa_dict[new_key] = str(record_dict.get(key).seq)

    return faa_dict


def genome_to_length_dict(genome_file):
    record_dict = SeqIO.index(genome_file, "fasta")

    length_dict = {}

    for key in record_dict.keys():

        length_dict[key] = len(record_dict.get(key).seq)

    return length_dict


# 最终文件格式要求参考 https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_external_gene_calls.txt
def transform(gff_file, faa_file, genome_file, tsv_file):


    # 去除评论行
    df = pd.read_table(gff_file, header=None, comment="#")
    df.columns = ["contig", "Source", "Feature", "start", "stop", "Score", "Strand", "Frame", "Attributes"]
    # 取feature为 gene或者pseudogene 的行
    gene_df = df[df["Feature"].isin(["gene", "pseudogene"])][["contig", "start", "stop", "Strand", "Attributes"]]

    gene_df["locus_tag"] = gene_df["Attributes"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series)["locus_tag"]

    # 取gene索引的下一行
    product_df = df.iloc[gene_df.index + 1,]

    gene_df = gene_df.drop(columns="Attributes").set_index("locus_tag", drop=False)

    product_attribute_df = product_df["Attributes"].str.split(";").apply(
        lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series).drop(columns=["ID", "Parent", "Name", "protein_id", "transl_table"]).set_index("locus_tag",
                                                                                                  drop=True)


    gene_df["direction"] = gene_df["Strand"].apply(lambda x: "f" if x == "+" else "r")
    gene_df = gene_df.drop(columns="Strand")
    # 'partial' (whether it is a complete gene call, or a partial one; must be 1 for partial calls, and 0 for complete calls)
    gene_df["partial"] = product_attribute_df["pseudo"].apply(lambda x: 1 if x=="true" else 0)
    # 'call_type' (1 if it is coding, 2 if it is noncoding, or 3 if it is unknown (only gene calls with call_type = 1 will have amino acid sequences translated))
    gene_df["call_type"] = product_attribute_df["gbkey"].apply(lambda x: 1 if x=="CDS" else 2)
    gene_df["source"] = "pgap"
    gene_df["version"] = "2023-10-03.build7061"
    # print(gene_df)

    faa_dict = faa_to_dict(faa_file)

    gene_df["aa_sequence"] = gene_df["locus_tag"].apply(lambda x: faa_dict.get(x, np.nan))
    gene_df["gene_callers_id"] = gene_df["locus_tag"].apply(lambda x: int(x.split("_")[-1]))

    gene_df = gene_df.set_index("gene_callers_id", drop=True)

    gene_df = gene_df.drop(columns="locus_tag")


    # 更改pgap的夸断点stop位置问题，pgag的这种处理方式会是stop的数值大于整个contig的长度，会导致anvio报错，会出现在contig的最后一个基因
    length_dict = genome_to_length_dict(genome_file)

    contig_name_list = gene_df['contig'].unique()

    for contig_name in contig_name_list:
        target_row = gene_df[gene_df["contig"] == contig_name].iloc[-1,]
        stop = target_row["stop"]
        index_name = target_row.name
        contig_length = length_dict[contig_name]
        if stop > contig_length:
            gene_df.at[index_name, "stop"] = contig_length


    gene_df.to_csv(tsv_file, sep="\t")


if __name__ == '__main__':
    transform(snakemake.input["gff"], snakemake.input["faa"], snakemake.input["genome"], snakemake.output[0])