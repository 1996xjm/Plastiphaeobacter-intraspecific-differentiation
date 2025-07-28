import os
import sys
from os import path
# 添加包搜索路径
sys.path.append(path.join(path.dirname(__file__), ".."))
from log import Log
import pandas as pd
from Bio import SeqIO
from io import StringIO

def get_value_count(col):
    split_res = col.dropna().str.split("_").apply(lambda x:"_".join(x[:-1]))
    return split_res.value_counts()

def save_single_copy_aas(single_copy_gene_df):
    os.makedirs(single_copy_aas_dir, exist_ok=True)
    all_records = SeqIO.index(all_faa_file, "fasta")
    for cluster_id,row in single_copy_gene_df.iterrows():
        record_list = []

        for locus_tag in row:
            record_list.append(all_records.get(locus_tag))
        with open(path.join(single_copy_aas_dir,f"{cluster_id}.faa"), "w") as f:
            SeqIO.write(record_list,f,"fasta")

def save_all_cluster_aas(mcl_df):
    # print(mcl_df)
    os.makedirs(all_cluster_aas_dir, exist_ok=True)
    all_records = SeqIO.index(all_faa_file, "fasta")
    for cluster_id,row in mcl_df.iterrows():
        record_list = []

        for locus_tag in row:
            if not pd.isna(locus_tag):
                record_list.append(all_records.get(locus_tag))
        with open(path.join(all_cluster_aas_dir,f"{cluster_id}.faa"), "w") as f:
            SeqIO.write(record_list,f,"fasta")


def save_single_copy_fnas(single_copy_gene_df):
    os.makedirs(single_copy_fnas_dir, exist_ok=True)
    all_records = SeqIO.index(all_fna_file, "fasta")
    for cluster_id,row in single_copy_gene_df.iterrows():
        record_list = []

        for locus_tag in row:
            record_list.append(all_records.get(locus_tag))
        with open(path.join(single_copy_fnas_dir,f"{cluster_id}.fna"), "w") as f:
            SeqIO.write(record_list,f,"fasta")


def save_single_copy_genes_info(single_copy_gene_df):
    os.makedirs(single_copy_genes_annotation_dir, exist_ok=True)
    os.makedirs(single_copy_sample_aas_dir, exist_ok=True)
    reorder_df = single_copy_gene_df.apply(lambda row: pd.Series(sorted(row)), axis=1)
    sample_id_list = ["_".join(i.split("_")[:2]) for i in reorder_df.iloc[0].to_list()]
    reorder_df.columns = sample_id_list
    reorder_df.reset_index(level=0,inplace=True)
    all_records = SeqIO.index(all_faa_file, "fasta")
    for sample_id in sample_id_list:
        annotation_df = pd.read_table(path.join(merge_annotation_output, f"{sample_id}_annotation.tsv"), index_col=0)
        sample_df = reorder_df[["cluster_id", sample_id]].set_index(sample_id)
        new_annotation_df = sample_df.join(annotation_df, how="left")
        locus_tag_list = sorted(new_annotation_df.index.tolist(), key=lambda s: int(s.split('_')[-1]))
        new_annotation_df.loc[locus_tag_list, :].to_csv(
            path.join(single_copy_genes_annotation_dir, f"{sample_id}_single_copy_genes_annotation.tsv"),
            sep="\t")

        record_list = [all_records.get(locus_tag) for locus_tag in locus_tag_list]
        with open(path.join(single_copy_sample_aas_dir,f"{sample_id}.faa"), "w") as f:
            SeqIO.write(record_list,f,"fasta")


def save_genome_specific_aas_blast_info(sample, locus_tag_list):
    out_dir = path.join(genome_specific_genes_blastp_info_dir, sample)
    relative_unique_gene_dir = path.join(out_dir,"relative_unique_gene") # 相对特异基因，在其他样本存在同源性低的blast结果
    os.makedirs(relative_unique_gene_dir, exist_ok=True)
    absolute_unique_gene_list = []

    mean_stringIO = StringIO()
    mean_stringIO.write("gene_id\tmean_effective_average_coverage\tmean_percent_of_self_bitscore\n")
    for gene_id in locus_tag_list:
        target_blast_df = blastp_res_df[blastp_res_df["query_id"] == gene_id].copy()
        # 如果subject_id都是以样本id起始， 就说明是特异的，只存在旁系同源基因
        if all(target_blast_df['subject_id'].str.startswith(sample)):
            absolute_unique_gene_list.append(gene_id)
            continue


        self_bitscore = target_blast_df[target_blast_df["subject_id"] == gene_id]["bit_score"]

        target_blast_df["effective_average_coverage"] = target_blast_df.apply(
            lambda row: row["aln_length"] * 2 / (row["q_length"] + row["s_length"]) * row["perc_id"] / 100, axis=1)

        target_blast_df["percent_of_self_bitscore"] = target_blast_df.apply(
            lambda row: row["bit_score"]/self_bitscore, axis=1)

        target_blast_df.to_csv(path.join(relative_unique_gene_dir,f"{gene_id}_blast.tsv"), sep="\t", index=None)

        without_self_df = target_blast_df[~(target_blast_df['query_id'].str.startswith(sample) & target_blast_df['subject_id'].str.startswith(sample))] # 滤掉自己和旁系同源基因
        # print(without_self_df)
        mean_stringIO.write(f"{gene_id}\t{without_self_df['effective_average_coverage'].mean()}\t{without_self_df['percent_of_self_bitscore'].mean()}\n")
    # print(mean_stringIO.getvalue())
    with open(path.join(out_dir, f"{sample}_mean_EAC_PSB.tsv"), "w") as f:
        f.write(mean_stringIO.getvalue())

    Log.info(f"absolute unique gene in {sample}: {len(absolute_unique_gene_list)}")

    annotation_df = pd.read_table(path.join(merge_annotation_output, f"{sample}_annotation.tsv"), index_col=0)
    absolute_unique_gene_annotation_df = annotation_df.loc[absolute_unique_gene_list, :]
    absolute_unique_gene_annotation_df.to_csv(path.join(out_dir, f"{sample}_absolute_unique_gene_annotation.tsv"),
                                  sep="\t")

def save_genome_specific_aas(mcl_df, gene_cluster_df):
    os.makedirs(genome_specific_aas_dir, exist_ok=True)
    os.makedirs(genome_specific_genes_annotation_dir, exist_ok=True)
    all_records = SeqIO.index(all_faa_file, "fasta")

    for sample in gene_cluster_df.columns:
        record_list = []
        locus_tag_list = mcl_df[(gene_cluster_df[sample] > 0) & (gene_cluster_df.drop(columns=sample).sum(axis=1) == 0)].stack().dropna().tolist()
        # 根据数字部分进行升序排序
        locus_tag_list = sorted(locus_tag_list, key=lambda s:int(s.split('_')[-1]))
        for locus_tag in locus_tag_list:
            record_list.append(all_records.get(locus_tag))
        with open(path.join(genome_specific_aas_dir, f"{sample}.faa"), "w") as f:
            SeqIO.write(record_list, f, "fasta")

        Log.info(f"specific aas in {sample}: {len(record_list)}")

        annotation_df = pd.read_table(path.join(merge_annotation_output, f"{sample}_annotation.tsv"), index_col=0)
        specific_annotation_df = annotation_df.loc[locus_tag_list, :]
        specific_annotation_df.to_csv(path.join(genome_specific_genes_annotation_dir, f"{sample}_specific_annotation.tsv"), sep="\t")
        new_strIO = StringIO()


        save_genome_specific_aas_blast_info(sample, locus_tag_list)


        # 存储chiplot 基因簇截断数据
        new_strIO.write("break_id\tsource\tbreakStart\tbreakEnd\n")
        for i in range(specific_annotation_df.shape[0]-1):
            current_end = specific_annotation_df.iloc[i,2]
            next_start = specific_annotation_df.iloc[i+1,1]

            current_seqname = specific_annotation_df.iloc[i,0]
            next_seqame = specific_annotation_df.iloc[i+1,0]

            margin = 1000

            # print(current_end, next_start , current_seqname)

            distance = 10000 # 跨度阈值/bp
            if (next_start - current_end) > distance and current_seqname == next_seqame:
                new_strIO.write(f"b_{i}\t{current_seqname}\t{current_end+margin}\t{next_start-margin}\n")
        new_strIO.seek(0)
        with open(path.join(genome_specific_genes_annotation_dir, f"{sample}_gene_cluster_break_list.tsv"), "w") as f:
            f.write(new_strIO.read())


def add_cluster_id_to_merge_annotation(mcl_df, gene_cluster_df):
    df = mcl_df.copy()
    df['cluster_id'] = df.index
    df_melted = df.melt(id_vars='cluster_id', value_name='gene_id').drop(columns=['variable']).dropna(axis=0)
    df_final = df_melted.set_index('gene_id')
    print("####dsfds", df_final)

    for sample in gene_cluster_df.columns:
        annotation_df = pd.read_table(path.join(merge_annotation_output, f"{sample}_annotation.tsv"), index_col=0)
        df_final.join(annotation_df, how="right").to_csv(path.join(merge_annotation_output, f"{sample}_annotation_with_clusterID.tsv"), sep="\t")


    return df_final

def transform():
    mcl_df = pd.read_table(mcl_clusters_file, header=None)

    mcl_df.index = [f"GC_{i + 1}" for i in mcl_df.index]
    mcl_df.index.name = "cluster_id"

    gene_cluster_df = mcl_df.T.apply(get_value_count, axis=0).fillna(0).T

    Log.info(f"total gene count: {gene_cluster_df.sum().sum()}")
    gene_cluster_df.to_csv(gene_cluster_count, sep="\t")


    single_copy_gene_df = mcl_df[(gene_cluster_df == 1).all(axis=1)].dropna(axis=1)
    core_gene_df = mcl_df[(gene_cluster_df >= 1).all(axis=1)]
    Log.info(f"total gene cluster count: {mcl_df.shape[0]}")
    Log.info(f"total core gene cluster count: {core_gene_df.shape[0]}")
    pair_split = snakemake.wildcards.sample.split("_vs_")
    pd.DataFrame({"pair1":[pair_split[0]], "pair2":[pair_split[1]], "value":[core_gene_df.shape[0] / mcl_df.shape[0] * 100]}).to_csv(core_genome_proportion,index=None, sep="\t")
    Log.info(f"total core gene count: {gene_cluster_df[(gene_cluster_df >= 1).all(axis=1)].sum().sum()}")

    # core_gene_df = pd.melt(core_gene_df, var_name='Variable', value_name='gene_id').dropna(axis=0).drop(columns="Variable")
    # core_gene_df["gene_type"] = "core gene"
    # core_gene_df = core_gene_df.set_index("gene_id")
    #
    # dispensable_gene_df = mcl_df[~(gene_cluster_df >= 1).all(axis=1)]
    # dispensable_gene_df = pd.melt(dispensable_gene_df, var_name='Variable', value_name='gene_id').dropna(axis=0).drop(
    #     columns="Variable")
    # dispensable_gene_df["gene_type"] = "accessory gene"
    # dispensable_gene_df = dispensable_gene_df.set_index("gene_id")
    #
    # pd.concat([core_gene_df, dispensable_gene_df], axis=0).to_csv(gene_type_table, sep="\t")
    #
    #
    #
    #
    # Log.info(f"total core single copy gene cluster count: {single_copy_gene_df.shape[0]}")
    # Log.info(f"average total single copy gene ratio in genome: {gene_cluster_df.sum().apply(lambda x:single_copy_gene_df.shape[0]/x).mean()}")
    #
    #
    # save_single_copy_aas(single_copy_gene_df)
    # save_all_cluster_aas(mcl_df)
    # save_single_copy_fnas(single_copy_gene_df)
    # save_single_copy_genes_info(single_copy_gene_df)

    # save_genome_specific_aas(mcl_df, gene_cluster_df)

    # gene_id_GC_id_df = add_cluster_id_to_merge_annotation(mcl_df, gene_cluster_df)
    #
    # pd.concat([pd.concat([core_gene_df, dispensable_gene_df], axis=0), gene_id_GC_id_df], axis=1).to_csv(gene_type_table, sep="\t")






if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        sys.stdout = f
        sys.stderr = f
        mcl_clusters_file = snakemake.input.mcl
        all_faa_file = snakemake.input.faa
        # all_fna_file = snakemake.input.fna
        blastp_res_file = snakemake.input.blastp_res
        gene_cluster_count = snakemake.output.gene_cluster_count
        core_genome_proportion = snakemake.output.core_genome_proportion
        # gene_type_table = snakemake.output.gene_type_table
        # single_copy_aas_dir = snakemake.output.single_copy_aas_dir
        # single_copy_sample_aas_dir = snakemake.output.single_copy_sample_aas_dir
        # all_cluster_aas_dir = snakemake.output.all_cluster_aas_dir
        # single_copy_fnas_dir = snakemake.output.single_copy_fnas_dir
        # single_copy_genes_annotation_dir = snakemake.output.single_copy_genes_annotation_dir
        # genome_specific_aas_dir = snakemake.output.genome_specific_aas_dir
        # genome_specific_genes_annotation_dir = snakemake.output.genome_specific_genes_annotation_dir
        # genome_specific_genes_blastp_info_dir = snakemake.output.genome_specific_genes_blastp_info_dir
        # merge_annotation_output = snakemake.params.merge_annotation_output

        blastp_res_df = pd.read_table(blastp_res_file, names=["query_id", "subject_id", "perc_id", "q_covs", "q_length", "s_length", "aln_length", "mismatches", "gaps", "e_val", "bit_score"])



        transform()