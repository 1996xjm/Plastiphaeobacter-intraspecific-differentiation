import os

import pandas as pd
from os import path




def find_gene():




    track_num = 0
    orthologous_gene_track_template_file_f = open(orthologous_gene_track_template_file, "w")
    for target_sample in all_smaples:
        if target_sample == ref_sample:
            continue
        final_ref_sample_gene_ids_list = []
        id_mapping_dict = {}
        with open(mcl_cluster_file, "r") as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                gene_ids_list = line.split("\t")
                ref_sample_gene_ids_list = [i for i in gene_ids_list if i.startswith(ref_sample)]
                target_sample_gene_ids_list = [i for i in gene_ids_list if i.startswith(target_sample)]
                if len(ref_sample_gene_ids_list) > 0 and len(target_sample_gene_ids_list) > 0:
                    final_ref_sample_gene_ids_list.extend(ref_sample_gene_ids_list)
                    for id in ref_sample_gene_ids_list:
                        id_mapping_dict[id] = (f"GC_{line_number}",target_sample_gene_ids_list[0])
        # print(len(final_ref_sample_gene_ids_list))

        ref_annotation_df = pd.read_table(ref_annotation, index_col=0)

        all_gene_df = ref_annotation_df[ref_annotation_df["Seqname"] == location]

        orthologous_gene_df = ref_annotation_df.loc[final_ref_sample_gene_ids_list, :]
        orthologous_gene_df = orthologous_gene_df[
            orthologous_gene_df["Seqname"] == location]

        not_orthologous_gene_df = all_gene_df[~all_gene_df.index.isin(orthologous_gene_df.index)]

        print(f"target_sample: {target_sample} total: {all_gene_df.shape[0]} orthologous_gene: {orthologous_gene_df.shape[0]} not_orthologous_gene: {not_orthologous_gene_df.shape[0]}")

        orthologous_gene_gff_file = f"{target_sample}_orthologous_genes.gff"
        orthologous_gene_gff_f = open(path.join(gff_out_dir, orthologous_gene_gff_file), "w")
        for i, row in orthologous_gene_df.iterrows():
            leftend = row["Start"]
            rightend = row["End"]
            print("\t".join([f"{location}", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                             f"Name={id_mapping_dict[i][0]}:{i}>{id_mapping_dict[i][1]}"]), file=orthologous_gene_gff_f)


        orthologous_gene_gff_f.close()

        not_orthologous_gene_gff_file = f"{target_sample}_not_orthologous_genes.gff"

        not_orthologous_gene_gff_f = open(path.join(gff_out_dir, not_orthologous_gene_gff_file), "w")

        for i, row in not_orthologous_gene_df.iterrows():
            leftend = row["Start"]
            rightend = row["End"]
            print("\t".join([f"{location}", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                             f"Name={i}"]), file=not_orthologous_gene_gff_f)

        not_orthologous_gene_gff_f.close()

        track_size = 10.0

        print("\t".join([f"{track_start_position - track_num * 0.05}", f"{track_size}" , "true", "true", "false", "false", "CDS", "null", "", "148:103:189", orthologous_gene_gff_file, path.abspath(gff_out_dir)]), file=orthologous_gene_track_template_file_f)
        print("\t".join([f"{track_start_position - track_num * 0.05}", f"{track_size}" , "true", "true", "false", "false", "CDS", "null", "", "204:204:204", not_orthologous_gene_gff_file, path.abspath(gff_out_dir)]), file=orthologous_gene_track_template_file_f)
        track_num += 1

    orthologous_gene_track_template_file_f.close()


if __name__ == '__main__':
    mcl_cluster_file = snakemake.input["mcl_cluster_file"]
    ref_annotation = snakemake.input["ref_annotation_file"]
    gff_out_dir = snakemake.params["gff_out_dir"]
    os.makedirs(gff_out_dir,exist_ok=True)
    track_start_position = snakemake.params["track_start_position"]
    location = snakemake.params["location"]
    orthologous_gene_track_template_file = snakemake.output[0]
    ref_sample = snakemake.wildcards["sample"]
    all_smaples = snakemake.params["all_smaples"]

    find_gene()