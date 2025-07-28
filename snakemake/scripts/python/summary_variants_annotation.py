import math
import os
import re
import sys
import pandas as pd
from log import Log
from os import path
from io import StringIO
from sklearn.decomposition import PCA

# https://bredagenetics.com/amino-acid-codes-symbols/
AA_property_tsv_str = """Amino acid	Code_1	Code_3	Physico-chemical property
Alanine	A	Ala	Hydrophobic side chain
lsoleucine	I	Ile	Hydrophobic side chain
Leucine	L	Leu	Hydrophobic side chain
Methionine	M	Met	Hydrophobic side chain
Phenylalanine	F	Phe	Hydrophobic side chain
Tryptophan	W	Trp	Hydrophobic side chain
Tyrosine	Y	Tyr	Hydrophobic side chain
Valine	V	Val	Hydrophobic side chain
Aspartic Acid	D	Asp	Negative charge
Glutamic Acid	E	Glu	Negative charge
Asparagine	N	Asn	Polar uncharged side chain
Glutamine	Q	Gln	Polar uncharged side chain
Serine	S	Ser	Polar uncharged side chain
Threonine	T	Thr	Polar uncharged side chain
Arginine	R	Arg	Positive charge
Histidine	H	His	Positive charge
Lysine	K	Lys	Positive charge
Cysteine	C	Cys	Special cases
Glicine	G	Gly	Special cases
Proline	P	Pro	Special cases"""

def base_changes(single_copy_gene_snp_ann_df, sample_name):


    # print(single_copy_gene_snp_ann_df)

    single_copy_gene_snp_ann_df["Base_changes"] = single_copy_gene_snp_ann_df["HGVS.c"].str.slice(-3)

    base_changes_df = single_copy_gene_snp_ann_df["TYPE"].groupby(
        by=[single_copy_gene_snp_ann_df["Gene_ID"],
            single_copy_gene_snp_ann_df["Base_changes"]]).count().unstack().fillna(
        0)

    all_genes_base_changes_df = pd.DataFrame(base_changes_df.sum()).T

    Log.info(f"base changes sum of {sample_name}: {base_changes_df.sum()}")

    all_genes_base_changes_df.columns = pd.MultiIndex.from_tuples(
        [tuple(c.split('>')) for c in all_genes_base_changes_df.columns])

    order = ["A", "G", "T", "C"]

    # 样本内所有基因的碱基替换矩阵
    all_genes_base_changes_4_4_matrix_df = all_genes_base_changes_df.stack(0).T[0].T.reindex(index=order,
                                                                                             columns=order).fillna(0)

    all_genes_base_changes_4_4_matrix_df.to_csv(
        path.join(base_changes_info_dir, f"{sample_name}_all_genes_base_changes_count_matrix.tsv"), sep="\t")

    (all_genes_base_changes_4_4_matrix_df / all_genes_base_changes_4_4_matrix_df.sum().sum() * 100).to_csv(
        path.join(base_changes_info_dir, f"{sample_name}_all_genes_base_changes_percent_matrix.tsv"), sep="\t")

    # creating a Multi-index df with two Rows index (Original Bases and Changed Bases)
    base_changes_df.columns = pd.MultiIndex.from_tuples([tuple(c.split('>')) for c in base_changes_df.columns])

    individual_gene_base_changes_4_4_matrix_df = base_changes_df.stack(0, future_stack=True)
    # 样本内单个基因的碱基替换矩阵
    individual_gene_base_changes_4_4_matrix_df = individual_gene_base_changes_4_4_matrix_df.reindex(
        columns=["A", "C", "G", "T"])


def SNV_codon_position(single_copy_gene_snp_ann_df, sample_name):

    single_copy_gene_snp_ann_df["SNV_codon_pos"] = single_copy_gene_snp_ann_df["CDS.pos_CDS.length"].str.split("/", expand=True)[0].astype(int).apply(lambda pos: pos%3 if pos%3 != 0 else 3)

    SNV_codon_position_df = single_copy_gene_snp_ann_df["TYPE"].groupby(by=[single_copy_gene_snp_ann_df["Gene_ID"],
                                                             single_copy_gene_snp_ann_df[
                                                                 "SNV_codon_pos"]]).count().unstack().fillna(0)
    SNV_codon_position_df.index.name = sample_name
    SNV_codon_position_df.to_csv(path.join(SNV_codon_position_info_dir,f"{sample_name}_SNV_codon_position.tsv"), sep="\t")
    SNV_codon_position_df.sum().to_csv(path.join(SNV_codon_position_info_dir,f"{sample_name}_SNV_codon_position_sum.tsv"), sep="\t")
    return SNV_codon_position_df


def SNV_variant_type(single_copy_gene_snp_ann_df, sample_name):
    """
    ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9417173/
    *387* Variants annotated as 'synonymous_variant', 'stop_retained_variant',
    *388* 'splice_region_variant&stop_retained_variant' or 'splice_region_variant&synonymous_variant'
    *389* by SnpEff were considered as the ‘synonymous’ variants in the same formula.
    :param single_copy_gene_snp_ann_df:
    :param sample_name:
    :return:
    """
    SNV_variant_type_df = single_copy_gene_snp_ann_df["TYPE"].groupby(by=[single_copy_gene_snp_ann_df["Gene_ID"],
                                                             single_copy_gene_snp_ann_df[
                                                                 "Annotation"]]).count().unstack()
    SNV_variant_type_df.index.name = sample_name
    SNV_variant_type_df[["missense_variant","synonymous_variant"]].to_csv(path.join(SNV_variant_missense_synonymous_info_dir,f"{sample_name}_SNV_variant_type_with_NA.tsv"), sep="\t")
    SNV_variant_type_df = SNV_variant_type_df.fillna(0)
    SNV_variant_type_df[["missense_variant","synonymous_variant"]].to_csv(path.join(SNV_variant_missense_synonymous_info_dir,f"{sample_name}_SNV_variant_type_no_NA.tsv"), sep="\t")
    SNV_variant_type_df.to_csv(path.join(SNV_variant_all_types_info_dir,f"{sample_name}_SNV_variant_type.tsv"), sep="\t")
    snv_variant_type_df_sum = SNV_variant_type_df.sum()
    snv_variant_type_df_sum.to_csv(path.join(SNV_variant_all_types_info_dir, f"{sample_name}_SNV_variant_type_sum.tsv"), sep="\t")
    snv_variant_type_df_sum.name = sample_name

    missense_variant_S = SNV_variant_type_df["missense_variant"]
    missense_variant_S.name = sample_name
    return snv_variant_type_df_sum,missense_variant_S


def AA_changes(single_copy_gene_snp_ann_df, AA_property_df, sample_name):
    # 只取错义突变的SNV
    missense_variant_df = single_copy_gene_snp_ann_df[single_copy_gene_snp_ann_df["Annotation"] == "missense_variant"].copy()
    missense_variant_df["AA_changes"] = missense_variant_df["HGVS.p"].apply(lambda x:re.sub(r'\d+', '>', x[2:]))
    print(missense_variant_df)



    AA_changes_df = missense_variant_df["TYPE"].groupby(
        by=[missense_variant_df["Gene_ID"],
            missense_variant_df["AA_changes"]]).count().unstack().fillna(
        0)

    AA_changes_df_sum = AA_changes_df.sum()

    all_genes_AA_changes_df = pd.DataFrame(AA_changes_df_sum).T

    Log.info(f'AA changes sum of {sample_name}: {all_genes_AA_changes_df.sum().sort_values()}')


    all_genes_AA_changes_df.columns = pd.MultiIndex.from_tuples(
        [tuple(c.split('>')) for c in all_genes_AA_changes_df.columns])

    all_genes_AA_changes_matrix_df = all_genes_AA_changes_df.stack(0).T[0].T.reindex(index=AA_property_df.index, columns=AA_property_df.index).fillna(0)
    all_genes_AA_changes_matrix_df.to_csv(path.join(AA_changes_info_dir,f"{sample_name}_all_genes_AA_changes.tsv"), sep="\t")

    return AA_changes_df_sum / AA_changes_df_sum.sum() * 100

# 氨基酸突变偏好
def AA_bias(single_copy_gene_snp_ann_df, AA_property_df, sample_name):
    # 只取错义突变的SNV
    missense_variant_df = single_copy_gene_snp_ann_df[single_copy_gene_snp_ann_df["Annotation"] == "missense_variant"].copy()
    missense_variant_df["AA_bias"] = missense_variant_df["HGVS.p"].apply(lambda x:x[-3:])




    AA_bias_df = missense_variant_df["TYPE"].groupby(
        by=[missense_variant_df["Gene_ID"],
            missense_variant_df["AA_bias"]]).count().unstack().fillna(
        0)

    # print(AA_bias_df)

    AA_bias_df_sum = AA_bias_df.sum()

    # print(AA_bias_df_sum)

    return AA_bias_df_sum / AA_bias_df_sum.sum() * 100


def get_loading_scores(count_df):
    # 实例化PCA  n_components 为最后保留的主成分个数 我们要画二维图 所以保留前两个主成分
    pca = PCA(n_components=2)
    # 数据预处理
    # print(df.values)
    pca.fit(count_df.values)
    # 获取降维后的数据
    X_dr = pca.transform(count_df.values)

    loading_scores = pd.DataFrame(pca.components_, index=["PC1", "PC2"], columns=count_df.columns).T
    loading_scores.index.name = "feature_id"
    loading_scores["point_distance"] = loading_scores.apply(lambda row: math.sqrt(row["PC1"] ** 2 + row["PC2"] ** 2),
                                                            axis=1)
    loading_scores = loading_scores.sort_values(by="point_distance", ascending=False)

    return loading_scores
def SNV_count_PCA(SNV_count_df, single_copy_gene_info_df):

    loading_scores = get_loading_scores(SNV_count_df)
    pd.concat([SNV_count_df.T,loading_scores],axis=1).sort_values(by="point_distance", ascending=False).to_csv(path.join(SNV_count_PCA_dir, "SNA_count_PCA_loading_scores_distance.tsv"), sep="\t")
    pd.concat([loading_scores, single_copy_gene_info_df], axis=1).to_csv(path.join(SNV_count_PCA_dir,"SNA_count_PCA_loading_scores_distance_with_function_annotation.tsv"),sep="\t")

    # print(loading_scores.sort_values(by="point_distance",ascending=False).to_csv(path.join(SNV_count_PCA_dir,"SNA_count_PCA_loading_scores_distance_with_function_annotation.tsv"),sep="\t"))
    return loading_scores


def summary():


    AA_property_df = pd.read_table(StringIO(AA_property_tsv_str), index_col=2)

    print(AA_property_df)

    single_copy_gene_info_df = pd.read_table(single_copy_gene_info_file, index_col=0)
    single_copy_gene_list = single_copy_gene_info_df.index.tolist()

    SNV_count_Series_list = []
    SNV_density_Series_list = []
    SNV_variant_type_Series_list = []
    missense_variant_Series_list = []
    AA_changes_Series_list = []
    AA_bias_Series_list = []
    SNV_codon_position_df_list = []
    for ann_vcf_file in ann_vcf_file_list:
        sample_name = path.basename(ann_vcf_file).replace("_ann.vcf", "")
        ann_vcf_df = pd.read_table(ann_vcf_file, comment="#", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown"])
        info_df = ann_vcf_df["INFO"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(pd.Series)
        ann_df = info_df["ANN"].str.split("|", expand=True).iloc[:,:14]
        ann_df.columns = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length", "AA.pos_AA.length"]
        ann_df["TYPE"] = info_df["TYPE"]
        # print(ann_df)



        single_copy_gene_ann_df = ann_df[ann_df["Gene_ID"].isin(single_copy_gene_list)].copy()
        single_copy_gene_snp_ann_df = ann_df[ann_df["Gene_ID"].isin(single_copy_gene_list) & (
                single_copy_gene_ann_df["TYPE"] == "snp")].copy()

        base_changes(single_copy_gene_snp_ann_df, sample_name)

        SNV_codon_position_df = SNV_codon_position(single_copy_gene_snp_ann_df, sample_name)
        SNV_codon_position_df_list.append(SNV_codon_position_df)

        snv_variant_type_df_sum,missense_variant_S = SNV_variant_type(single_copy_gene_snp_ann_df, sample_name)

        SNV_variant_type_Series_list.append(snv_variant_type_df_sum)
        missense_variant_Series_list.append(missense_variant_S)



        AA_changes_percentage = AA_changes(single_copy_gene_snp_ann_df, AA_property_df, sample_name)
        AA_changes_percentage.name = sample_name
        AA_changes_Series_list.append(AA_changes_percentage)

        AA_bias_percentage = AA_bias(single_copy_gene_snp_ann_df, AA_property_df, sample_name)
        AA_bias_percentage.name = sample_name

        AA_bias_Series_list.append(AA_bias_percentage)






        Log.info(f'mutation type of {sample_name}:')
        print(single_copy_gene_ann_df["TYPE"].value_counts().to_csv(sep="\t"))
        # 计算每个gene ID突变类型的数目
        mutation_type_df = single_copy_gene_ann_df["Allele"].groupby(
            by=[single_copy_gene_ann_df["Gene_ID"], single_copy_gene_ann_df["TYPE"]]).count().unstack().fillna(0)

        mutation_type_df["gene_length"] = single_copy_gene_info_df["End"] - single_copy_gene_info_df["Start"] + 1
        mutation_type_df["SNV_density"] = mutation_type_df["snp"]/mutation_type_df["gene_length"] * 100
        # print(mutation_type_df)
        SNV_count_Series = mutation_type_df["snp"]
        SNV_count_Series.name = sample_name
        SNV_count_Series_list.append(SNV_count_Series)

        SNV_density_Series = mutation_type_df["SNV_density"]
        SNV_density_Series.name = sample_name
        SNV_density_Series_list.append(SNV_density_Series)

        # print(mutation_type_df)
    SNV_count_df = pd.concat(SNV_count_Series_list, axis=1)
    Log.info(f"total single copy core genes that have SNVs: {SNV_count_df.shape[0]}")
    SNV_count_df.to_csv(path.join(SNV_info_dir, "SNV_count_with_NA.tsv"), sep="\t")
    SNV_count_df.fillna(0).to_csv(path.join(SNV_info_dir, "SNV_count_no_NA.tsv"), sep="\t")
    SNV_count_df.fillna(0).T.to_csv(path.join(SNV_info_dir, "SNV_count_no_NA_PCA.tsv"), sep="\t")

    SNV_count_PCA_loadind_scores_df = SNV_count_PCA(SNV_count_df.fillna(0).T, single_copy_gene_info_df)



    Log.info("These gene IDs hava no SNVs across all samples:")

    print(single_copy_gene_info_df.index[~single_copy_gene_info_df.index.isin(SNV_count_df.index)].tolist())

    SNV_density_df = pd.concat(SNV_density_Series_list, axis=1)
    SNV_density_df.to_csv(path.join(SNV_info_dir, "SNV_density_with_NA.tsv"), sep="\t")
    SNV_density_df.fillna(0).to_csv(path.join(SNV_info_dir, "SNV_density_no_NA.tsv"), sep="\t")
    # print(SNV_density_df)

    pd.concat(SNV_variant_type_Series_list, axis=1).fillna(0).T.to_csv(path.join(SNV_variant_type_info_dir, "all_samples_SNV_variant_type_sum.tsv"), sep="\t")
    pd.concat(missense_variant_Series_list,axis=1).fillna(0).T.to_csv(path.join(SNV_variant_type_info_dir, "all_samples_missense_variant_count_PCA.tsv"), sep="\t")
    AA_changes_percentage_df = pd.concat(AA_changes_Series_list, axis=1).fillna(0)
    AA_changes_percentage_df = AA_changes_percentage_df.loc[AA_changes_percentage_df.mean(axis=1).sort_values(ascending=False).index]
    AA_changes_percentage_df.T.to_csv(path.join(AA_changes_info_dir, "AA_changes_percentage_in_each_sample.tsv"), sep="\t")
    AA_changes_percentage_df.T.iloc[:,:25].to_csv(path.join(AA_changes_info_dir, "AA_changes_percentage_in_each_sample_top25.tsv"), sep="\t")



    AA_bias_percentage_df = pd.concat(AA_bias_Series_list, axis=1).fillna(0)
    AA_bias_percentage_df = AA_bias_percentage_df.loc[
        AA_bias_percentage_df.mean(axis=1).sort_values(ascending=False).index].T
    AA_bias_percentage_df.to_csv(path.join(AA_changes_info_dir, "AA_bias_percentage_in_each_sample.tsv"), sep="\t")

    def get_pos_Series(SNV_codon_position_df_list,pos=1):
        res_list = []

        for df in SNV_codon_position_df_list:
            pos_Series = df[pos]
            pos_Series.name = df.index.name
            res_list.append(pos_Series)

        return res_list

    def get_pos_ratio_Series(SNV_codon_position_df_list,pos=1):
        res_list = []

        for df in SNV_codon_position_df_list:
            ratio_df = df.div(df.sum(axis=1), axis=0)
            pos_Series = ratio_df[pos]
            pos_Series.name = df.index.name
            res_list.append(pos_Series)

        return res_list




    gene_id_set_dict = {}
    gene_id_set_dict["all_pos"] = SNV_count_PCA_loadind_scores_df.index.tolist()[:150]
    for i in range(1,4):
        res_df = pd.concat(get_pos_Series(SNV_codon_position_df_list, i), axis=1).fillna(0).T
        res_df.to_csv(path.join(SNV_codon_position_info_dir, f"all_samples_SNV_codon_position_{i}_count.tsv"), sep="\t")
        gene_id_set_dict[f"pos_{i}"]= get_loading_scores(res_df).index.tolist()[:150]

    pd.DataFrame(gene_id_set_dict).to_csv(path.join(SNV_codon_position_info_dir, f"geneID_SNV_codon_position_loading_scores_top150.tsv"), sep="\t", index=False)


    for i in range(1,4):
        pd.concat(get_pos_ratio_Series(SNV_codon_position_df_list,i),axis=1).fillna(0).T.to_csv(path.join(SNV_codon_position_info_dir,f"all_samples_SNV_codon_position_{i}_ratio.tsv"), sep="\t")






    # print(ann_df.to_csv(path.join(SNV_count_dir, "test.tsv"), sep="\t"))

if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        # sys.stdout = f
        # sys.stderr = f
        ann_vcf_file_list = snakemake.input.ann_vcf_list
        single_copy_gene_info_file = snakemake.input.single_copy_gene_info

        SNV_info_dir = snakemake.output.SNV_info_dir
        SNV_count_PCA_dir = path.join(SNV_info_dir, "SNV_count_PCA")

        base_changes_info_dir = snakemake.output.base_changes_info_dir
        SNV_codon_position_info_dir = snakemake.output.SNV_codon_position_info_dir
        SNV_variant_type_info_dir = snakemake.output.SNV_variant_type_info_dir
        SNV_variant_all_types_info_dir = path.join(SNV_variant_type_info_dir,"all_types")
        SNV_variant_missense_synonymous_info_dir = path.join(SNV_variant_type_info_dir,"missense_synonymous_variant")

        AA_changes_info_dir = snakemake.output.AA_changes_info_dir
        os.makedirs(SNV_count_PCA_dir, exist_ok=True)
        os.makedirs(base_changes_info_dir, exist_ok=True)
        os.makedirs(SNV_codon_position_info_dir, exist_ok=True)

        os.makedirs(SNV_variant_all_types_info_dir, exist_ok=True)
        os.makedirs(SNV_variant_missense_synonymous_info_dir, exist_ok=True)

        os.makedirs(AA_changes_info_dir, exist_ok=True)

        summary()