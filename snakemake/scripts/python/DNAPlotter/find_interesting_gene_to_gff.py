import os

import pandas as pd
from os import path
import json




def find_transposase_gene():

    ref_annotation_df = pd.read_table(ref_annotation, index_col=0)
    transposase_df = ref_annotation_df[
        ref_annotation_df["product"].str.contains("transposase", case=False) |
        ref_annotation_df["bakta_Product"].str.contains("transposase", case=False) |
        ref_annotation_df["rast_Product"].str.contains("transposase", case=False)
    ]

    transposase_df = transposase_df[
        transposase_df["Seqname"] == location]

    # print(transposase_df)

    transposase_gene_gff_f = open(transposase_gene_gff_file, "w")

    for i, row in transposase_df.iterrows():
        leftend = row["Start"]
        rightend = row["End"]
        print("\t".join([f"{location}", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                         f"Name={i}"]), file=transposase_gene_gff_f)

    transposase_gene_gff_f.close()


def find_recombinase_integrase_gene():

    ref_annotation_df = pd.read_table(ref_annotation, index_col=0)
    integrase_df = ref_annotation_df[
        ref_annotation_df["product"].str.contains("integrase", case=False) |
        ref_annotation_df["product"].str.contains("recombinase", case=False) |
        ref_annotation_df["bakta_Product"].str.contains("integrase", case=False) |
        ref_annotation_df["bakta_Product"].str.contains("recombinase", case=False) |
        ref_annotation_df["rast_Product"].str.contains("integrase", case=False) |
        ref_annotation_df["rast_Product"].str.contains("recombinase", case=False)

        ]

    integrase_df = integrase_df[
        integrase_df["Seqname"] == location]

    # print(integrase_df[["product","bakta_Product","rast_Product"]].to_csv(sep="\t"))

    integrase_gene_gff_f = open(recombinase_integrase_gene_gff_file, "w")

    for i, row in integrase_df.iterrows():
        leftend = row["Start"]
        rightend = row["End"]
        print("\t".join([f"{location}", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                         f"Name={i}"]), file=integrase_gene_gff_f)

    integrase_gene_gff_f.close()



def find_phage_region():
    with open(chromosome_phage_regions, "r") as f:
        chromosome_GTA_gff_file_f = open(chromosome_GTA_gff_file, "w")
        phage_region_intact_gff_f = open(chromosome_phage_regions_intact_gff_file, "w")
        phage_region_questionable_gff_f = open(chromosome_phage_regions_questionable_gff_file, "w")
        phage_region_incomplete_gff_f = open(chromosome_phage_regions_incomplete_gff_file, "w")
        f_dict = {
            "GTA":chromosome_GTA_gff_file_f,
            "intact":phage_region_intact_gff_f,
            "questionable":phage_region_questionable_gff_f,
            "incomplete":phage_region_incomplete_gff_f
        }
        for region in json.load(f):
            leftend = region["start"]
            rightend = region["stop"]
            region_index = region["region"]
            region_completeness = region["completeness"]

            print("\t".join([f"{location}", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                             f"Name=PR_{region_index}: {region_completeness}"]), file=f_dict[region_completeness])

        chromosome_GTA_gff_file_f.close()
        phage_region_intact_gff_f.close()
        phage_region_questionable_gff_f.close()
        phage_region_incomplete_gff_f.close()




if __name__ == '__main__':
    seq_type = snakemake.params["seq_type"]
    location = snakemake.params["location"]

    ref_annotation = snakemake.input["ref_annotation_file"]
    transposase_gene_gff_file = snakemake.output["transposase_gene_gff_file"]
    recombinase_integrase_gene_gff_file = snakemake.output["recombinase_integrase_gene_gff_file"]


    if seq_type == "chromosome":
        chromosome_phage_regions = snakemake.input["chromosome_phage_regions"]
        chromosome_GTA_gff_file = snakemake.output["chromosome_GTA_gff_file"]
        chromosome_phage_regions_intact_gff_file = snakemake.output["chromosome_phage_regions_intact_gff_file"]
        chromosome_phage_regions_questionable_gff_file = snakemake.output["chromosome_phage_regions_questionable_gff_file"]
        chromosome_phage_regions_incomplete_gff_file = snakemake.output["chromosome_phage_regions_incomplete_gff_file"]
        find_phage_region()




    find_transposase_gene()
    find_recombinase_integrase_gene()
