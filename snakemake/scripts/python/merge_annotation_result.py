
import pandas as pd



def get_position_identifier(row):

    return f"{row['Seqname']}..{row['Start']}..{row['End']}..{row['Strand']}"

def get_position_identifier_rast(row):

    if row["Strand"] == "+":
        return f"{row['Contig']}..{row['Start']}..{row['End']}..{row['Strand']}"
    else:
        return f"{row['Contig']}..{row['End']}..{row['Start']}..{row['Strand']}"


def pgapGFF_bakta_to_tsv(gff_file, bakta_tsv, tsv_file):
    """
    用bakta自己预测基因并注释的结果
    :param gff_file:
    :param bakta_tsv:
    :param tsv_file:
    :return:
    """

    # 去除评论行
    df = pd.read_table(gff_file, header=None, comment="#")
    df.columns = ["Seqname", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]
    # 取feature为 gene或者pseudogene 的行
    gene_df = df[df["Feature"].isin(["gene", "pseudogene"])][["Seqname", "Start", "End", "Strand", "Attributes"]]

    gene_df["locus_tag"] = gene_df["Attributes"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series)["locus_tag"]

    gene_df["Position_identifier"] = gene_df.apply(get_position_identifier, axis=1)


    # 取gene索引的下一行
    product_df = df.iloc[gene_df.index + 1,]

    gene_df = gene_df.drop(columns="Attributes").set_index("locus_tag", drop=False)

    product_attribute_df = product_df["Attributes"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series).drop(columns=["ID", "Parent", "Name", "protein_id", "transl_table"]).set_index("locus_tag", drop=True)

    product_attribute_df["Note"] = product_attribute_df["Note"].str.replace("%3B",";")


    pgap_df = pd.concat([gene_df, product_attribute_df], axis=1)


    bakta_df = pd.read_table(bakta_tsv, header=None, comment="#")
    bakta_df.columns = ["Seqname", "Type", "Start", "End", "Strand", "bakta_Locus_Tag", "bakta_Gene", "bakta_Product", "bakta_DbXrefs"]

    bakta_df["Position_identifier"] = bakta_df.apply(get_position_identifier, axis=1)

    bakta_df = bakta_df.drop(columns=["Seqname", "Type", "Start", "End", "Strand"])


    # 不是所有的 Position_identifier都能匹配，因为某些基因（存在多个起始密码子选择）在两个软件预测出的位置不一样，bakta的策略是贪婪模式
    final_df = pd.merge(pgap_df, bakta_df, left_on="Position_identifier", right_on="Position_identifier", how="left")
    # print(final_df[["product","bakta_Product"]])

    final_df = final_df.set_index("locus_tag",drop=True)

    final_df.to_csv(tsv_file,sep="\t")



def pgapGFF_bakta_proteins_rast_to_tsv(gff_file, bakta_tsv, rast_tsv, tsv_file):
    """
    用pgap的蛋白序列给bakta_proteins注释的结果
    :param gff_file:
    :param bakta_tsv:
    :param tsv_file:
    :return:
    """

    # 去除评论行
    df = pd.read_table(gff_file, header=None, comment="#")
    df.columns = ["Seqname", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]
    # 取feature为 gene或者pseudogene 的行
    gene_df = df[df["Feature"].isin(["gene", "pseudogene"])][["Seqname", "Start", "End", "Strand", "Attributes"]]

    gene_df["locus_tag"] = gene_df["Attributes"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series)["locus_tag"]

    gene_df["Position_identifier"] = gene_df.apply(get_position_identifier, axis=1)



    # 取gene索引的下一行
    product_df = df.iloc[gene_df.index + 1,]

    gene_df = gene_df.drop(columns="Attributes").set_index("locus_tag", drop=False)

    product_attribute_df = product_df["Attributes"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(
        pd.Series).drop(columns=["ID", "Parent", "Name", "protein_id", "transl_table"]).set_index("locus_tag", drop=True)

    product_attribute_df["Note"] = product_attribute_df["Note"].str.replace("%3B",";")


    pgap_df = pd.concat([gene_df, product_attribute_df], axis=1)


    bakta_df = pd.read_table(bakta_tsv, comment="#")
    # print(bakta_df)


    bakta_df = bakta_df.drop(columns=["UniParc", "UniRef"])
    bakta_df.columns = ["ID", "bakta_Length", "bakta_Gene", "bakta_Product", "bakta_EC", "bakta_GO", "bakta_COG", "bakta_RefSeq"]
    bakta_df["locus_tag"] = bakta_df["ID"].str.split("|", expand=True)[2]
    bakta_df = bakta_df.set_index("locus_tag", drop=True)

    bakta_df = bakta_df.drop(columns=["ID"])

    final_df = pd.concat([pgap_df,bakta_df], axis=1)

    rast_df = pd.read_table(rast_tsv)

    rast_df["Position_identifier"] = rast_df.apply(get_position_identifier_rast, axis=1)

    rast_df = rast_df[["Loctag ID", "Function", "Subsystem", "Position_identifier"]]
    rast_df.columns = ["rast_locus_tag", "rast_Product", "rast_Subsystem", "Position_identifier"]

    final_df = pd.merge(final_df, rast_df, how="left", left_on="Position_identifier", right_on="Position_identifier")

    final_df = final_df.drop(columns="Position_identifier").set_index("locus_tag", drop=True)





    final_df.to_csv(tsv_file,sep="\t")


if __name__ == '__main__':
    pgapGFF_bakta_proteins_rast_to_tsv(snakemake.input["pgap_gff"], snakemake.input["bakta_tsv"], snakemake.input["rast_tsv"], snakemake.output[0])