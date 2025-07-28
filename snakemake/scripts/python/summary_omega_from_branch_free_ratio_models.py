import math
import os
from os import path
import re
import subprocess

import numpy as np
from sklearn.decomposition import PCA
from log import Log
import pandas as pd

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

# def summary():
#     Series_list = []
#     for file in codeml_out_list:
#         dNdS_Series_list = []
#         grep_list = [("dN", "dN tree:"), ("dS", "dS tree:")]
#         for k, v in grep_list:
#             output, error = subprocess.Popen(f"grep -A 1 '{v}' {file}", stdout=subprocess.PIPE, shell=True).communicate()
#             tree_str = output.decode().strip().split("\n")[1]
#
#             # print(tree_str)
#             # Regular expression to match leaf names and values
#             matches = re.findall(r'(\w+): (\d+\.\d+)', tree_str)
#
#             GC_dict = {}
#
#             for match in matches:
#                 GC_dict["_".join(match[0].split("_")[:2])] = match[1]
#                 # print(f"Leaf Name: {match[0]}, Value: {match[1]}")
#
#             S = pd.Series(GC_dict).astype(float)
#             S.name = k
#             dNdS_Series_list.append(S)
#
#
#         output, error = subprocess.Popen(f"grep -A 4 'dN & dS for each branch' {file}", stdout=subprocess.PIPE, shell=True).communicate()
#         dndn_branch_row = re.sub(r" +", "\t", output.decode().strip().split("\n")[-1]).split("\t")
#
#         N = float(dndn_branch_row[3])
#         S = float(dndn_branch_row[4])
#         # print(N,S)
#
#         # if the observed count of sites (pN or pS) less than 1, this branch will be regarded as extremely purifying selection
#         omega_Series = pd.concat(dNdS_Series_list, axis=1).apply(
#             lambda row: 0 if row["dN"] * N < 1 or row["dS"] * S < 1 else
#             row["dN"] / row["dS"], axis=1)
#
#         omega_Series.name = file.split("/")[-2]
#         Series_list.append(omega_Series)
#     df = pd.concat(Series_list, axis=1)
#     df.T.to_csv(omega_matrix, sep="\t")
#     # drop columns with the same value
#     PCA_df = df.loc[:, df.nunique() > 1]
#     PCA_df.to_csv(omega_PCA, sep="\t")
#
#     PCA_less_1 = PCA_df.loc[:, ~(PCA_df > 1).any()]
#     # PCA_less_1.to_csv(omega_less_1_PCA, sep="\t")
#     #
#     # PCA_less_1.loc[:, (PCA_less_1 == 0).sum() < PCA_less_1.shape[0] / 2].to_csv(omega_less_1_and_no_half_zero_PCA, sep="\t")
#     # PCA_less_1.loc[:, ~(PCA_less_1 == 0).any()].to_csv(omega_less_1_and_no_zero_PCA, sep="\t")
#
#     SNV_density_df = pd.read_table(SNV_density, index_col=0)
#
#     PCA_SNV_density_df = df.loc[:, SNV_density_df[SNV_density_df["SNV_density"] >= 0.1].index.tolist()]
#     # PCA_SNV_density_df = PCA_SNV_density_df.loc[:, PCA_SNV_density_df.nunique() > 1]
#     PCA_SNV_density_df.to_csv(omega_SNV_density_greater_10percent_PCA, sep="\t")
#     loading_scores = get_loading_scores(PCA_SNV_density_df)
#     print(loading_scores)
#     omega_SNV_density_greater_10percent_matrix_df = PCA_SNV_density_df.loc[:, loading_scores.index].T
#     omega_SNV_density_greater_10percent_matrix_df.to_csv(omega_SNV_density_greater_10percent_matrix, sep="\t")
#     print(pd.read_table(single_copy_gene_info_file).set_index("cluster id").loc[omega_SNV_density_greater_10percent_matrix_df.index,["bakta_Product","PT42_8"]].to_csv(sep="\t"))


def summary2():
    """
    use log10 to standardize the dN/dS value
    :return:
    """
    Series_list = []
    for file in codeml_out_list:
        output, error = subprocess.Popen(f"grep -A 1 'w ratios as node labels:' {file}", stdout=subprocess.PIPE, shell=True).communicate()
        tree_str = output.decode().strip().split("\n")[1]

        # print(tree_str)
        # Regular expression to match leaf names and values
        matches = re.findall(r'(\w+) #(\d+\.?\d+)', tree_str)

        GC_dict = {}

        for match in matches:
            GC_dict["_".join(match[0].split("_")[:2])] = match[1]
            # print(f"Leaf Name: {match[0]}, Value: {match[1]}")

        omega_Series = pd.Series(GC_dict).astype(float)
        omega_Series.name = file.split("/")[-2]
        Series_list.append(omega_Series)




    df = pd.concat(Series_list, axis=1)
    df_log10 = np.log10(df)

    df.T.to_csv(omega_matrix, sep="\t")

    # drop columns with the same value
    PCA_df = df_log10.loc[:, df_log10.nunique() > 1]
    PCA_df.to_csv(omega_PCA, sep="\t")

    S1 = set(get_loading_scores(PCA_df).index.tolist()[:265])
    print(S1)



    SNV_density_df = pd.read_table(SNV_density, index_col=0)
    df.loc[:, SNV_density_df[SNV_density_df["SNV_density"] >= 0.1].index.tolist()].T.to_csv(omega_SNV_density_greater_10percent_matrix, sep="\t")
    # SNV密度太低会导致计算出来的值要么是负无穷，要么是正无穷
    PCA_SNV_density_df = df_log10.loc[:, SNV_density_df[SNV_density_df["SNV_density"] >= 0.1].index.tolist()]
    # PCA_SNV_density_df = PCA_SNV_density_df.loc[:, PCA_SNV_density_df.nunique() > 1]
    PCA_SNV_density_df.to_csv(omega_SNV_density_greater_10percent_PCA, sep="\t")
    loading_scores = get_loading_scores(PCA_SNV_density_df)
    print(S1.intersection(set(loading_scores.index.tolist())))
    omega_SNV_density_greater_10percent_matrix_df = PCA_SNV_density_df.loc[:, loading_scores.index].T
    omega_SNV_density_greater_10percent_matrix_df.to_csv(omega_SNV_density_greater_10percent_matrix_log10, sep="\t")
    # print(pd.read_table(single_copy_gene_info_file).set_index("cluster id").loc[omega_SNV_density_greater_10percent_matrix_df.index,["bakta_Product","PT42_8"]].to_csv(sep="\t"))
    # (0.05, 0.1)
    PCA_SNV_density_2_df = df_log10.loc[:, SNV_density_df[(SNV_density_df["SNV_density"] >= 0.05) & (
                SNV_density_df["SNV_density"] < 0.1)].index.tolist()]
    print(get_loading_scores(PCA_SNV_density_2_df).to_csv(sep="\t"))




if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        # sys.stdout = f
        # sys.stderr = f


        codeml_out_list = snakemake.input.codeml_out_list
        SNV_density = snakemake.input.SNV_density
        single_copy_gene_info_file = snakemake.input.single_copy_gene_info_file

        omega_matrix = snakemake.output.omega_matrix
        omega_PCA = snakemake.output.omega_PCA

        omega_SNV_density_greater_10percent_PCA = snakemake.output.omega_SNV_density_greater_10percent_PCA
        omega_SNV_density_greater_10percent_matrix_log10 = snakemake.output.omega_SNV_density_greater_10percent_matrix_log10
        omega_SNV_density_greater_10percent_matrix = snakemake.output.omega_SNV_density_greater_10percent_matrix




        summary2()