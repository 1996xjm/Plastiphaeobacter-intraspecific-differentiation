import math
import os
import re
import sys
import pandas as pd
from log import Log
from os import path
from io import StringIO


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

def AA_changes(single_copy_gene_snp_ann_df, AA_property_df, sample_name):
    # 只取错义突变的SNV
    missense_variant_df = single_copy_gene_snp_ann_df[single_copy_gene_snp_ann_df["Annotation"] == "missense_variant"].copy()
    missense_variant_df["AA_changes"] = missense_variant_df["HGVS.p"].apply(lambda x:re.sub(r'\d+', '>', x[2:]))
    # print(missense_variant_df)



    AA_changes_df = missense_variant_df["TYPE"].groupby(
        by=[missense_variant_df["Gene_ID"],
            missense_variant_df["AA_changes"]]).count().unstack().fillna(
        0)

    print(AA_changes_df)

    AA_changes_df_sum = AA_changes_df.sum()

    # all_genes_AA_changes_df = pd.DataFrame(AA_changes_df_sum).T
    #
    # Log.info(f'AA changes sum of {sample_name}: {all_genes_AA_changes_df.sum().sort_values()}')
    #
    #
    # all_genes_AA_changes_df.columns = pd.MultiIndex.from_tuples(
    #     [tuple(c.split('>')) for c in all_genes_AA_changes_df.columns])
    #
    # all_genes_AA_changes_matrix_df = all_genes_AA_changes_df.stack(0).T[0].T.reindex(index=AA_property_df.index, columns=AA_property_df.index).fillna(0)
    # all_genes_AA_changes_matrix_df.to_csv(path.join(AA_changes_info_dir,f"{sample_name}_all_genes_AA_changes.tsv"), sep="\t")

    return AA_changes_df_sum



def summary():


    AA_property_df = pd.read_table(StringIO(AA_property_tsv_str), index_col=2)

    # print(AA_property_df)



    SNV_count_Series_list = []
    SNV_density_Series_list = []
    SNV_variant_type_Series_list = []
    missense_variant_Series_list = []
    AA_changes_Series_list = []
    AA_bias_Series_list = []
    SNV_codon_position_df_list = []
    for ann_vcf_file in ann_vcf_file_list:
        sample_name = path.basename(ann_vcf_file).replace("_ann.vcf", "")
        single_copy_gene_info_df = pd.read_table(path.join(single_copy_gene_ann_dir,f"{sample_name}_single_copy_genes_annotation.tsv"), index_col=0)
        single_copy_gene_list = single_copy_gene_info_df.index.tolist()
        ann_vcf_df = pd.read_table(ann_vcf_file, comment="#", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown"])
        info_df = ann_vcf_df["INFO"].str.split(";").apply(lambda x: dict(item.split("=") for item in x)).apply(pd.Series)
        ann_df = info_df["ANN"].str.split("|", expand=True).iloc[:,:14]
        ann_df.columns = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length", "AA.pos_AA.length"]
        ann_df["TYPE"] = info_df["TYPE"]
        # print(ann_df)



        single_copy_gene_ann_df = ann_df[ann_df["Gene_ID"].isin(single_copy_gene_list)].copy()
        single_copy_gene_snp_ann_df = ann_df[ann_df["Gene_ID"].isin(single_copy_gene_list) & (
                single_copy_gene_ann_df["TYPE"] == "snp")].copy()





        AA_changes_df = AA_changes(single_copy_gene_snp_ann_df, AA_property_df, sample_name)
        AA_changes_df.name = sample_name
        AA_changes_Series_list.append(AA_changes_df)
        #
        # AA_bias_percentage = AA_bias(single_copy_gene_snp_ann_df, AA_property_df, sample_name)
        # AA_bias_percentage.name = sample_name
        #
        # AA_bias_Series_list.append(AA_bias_percentage)






        Log.info(f'mutation type of {sample_name}:')
        print(single_copy_gene_ann_df["TYPE"].value_counts().to_csv(sep="\t"))
        # 计算每个gene ID突变类型的数目
        mutation_type_df = single_copy_gene_ann_df["Allele"].groupby(
            by=[single_copy_gene_ann_df["Gene_ID"], single_copy_gene_ann_df["TYPE"]]).count().unstack().fillna(0)

        mutation_type_df["gene_length"] = single_copy_gene_info_df["End"] - single_copy_gene_info_df["Start"] + 1
        mutation_type_df["SNV_density (Kbp)"] = mutation_type_df["snp"]/mutation_type_df["gene_length"] * 1000
        mutation_type_df["cluster_id"] = single_copy_gene_info_df["cluster_id"]
        mutation_type_df = mutation_type_df.set_index("cluster_id")
        # print(mutation_type_df)
        SNV_count_Series = mutation_type_df["snp"]
        SNV_count_Series.name = sample_name
        SNV_count_Series_list.append(SNV_count_Series)

        SNV_density_Series = mutation_type_df["SNV_density (Kbp)"]
        SNV_density_Series.name = sample_name
        SNV_density_Series_list.append(SNV_density_Series)

    SNV_count_df = pd.concat(SNV_count_Series_list, axis=1).fillna(0)
    SNV_count_df["Average count"] = SNV_count_df.mean(axis=1)
    SNV_count_df.to_csv(snakemake.output.SNV_count_file, sep="\t")

    SNV_density_df = pd.concat(SNV_density_Series_list, axis=1).fillna(0)
    SNV_density_df["Average density"] = SNV_density_df.mean(axis=1)
    SNV_density_df.to_csv(snakemake.output.SNV_density_file, sep="\t")

    AA_changes_percentage_df = pd.concat(AA_changes_Series_list, axis=1).fillna(0)
    AA_changes_percentage_df["Percentage (%)"] = AA_changes_percentage_df.sum(axis=1)/AA_changes_percentage_df.sum().sum() * 100
    AA_changes_percentage_df.to_csv(snakemake.output.AA_changes_percentage_file,sep="\t")


if __name__ == '__main__':
    ann_vcf_file_list = snakemake.input.ann_vcf_list
    single_copy_gene_ann_dir = snakemake.params.single_copy_gene_ann_dir
    # print(ann_vcf_file_list)
    summary()