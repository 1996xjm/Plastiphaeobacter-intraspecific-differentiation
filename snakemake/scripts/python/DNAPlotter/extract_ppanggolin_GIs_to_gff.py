from os import path
import pandas as pd
import re

def extract_GIs():
    GIs_df = pd.read_table(ppanggolin_GIs_file,index_col=0)
    GI_index = 1
    GI_gff_f = open(gff_out_path, "w")
    for i, row in GIs_df[GIs_df["genome"]==ref_sample].sort_values(by="start").iterrows():
        leftend = row["start"]
        rightend = row["stop"]
        print("\t".join([f"{ref_sample}_chromosome_1", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".",
                         f"Name={ref_sample}_GI_{GI_index}"]), file=GI_gff_f)
        GI_index += 1

if __name__ == '__main__':

    ref_sample = snakemake.wildcards["sample"]

    ppanggolin_GIs_file = snakemake.input["ppanggolin_GIs_file"]
    gff_out_path = snakemake.output[0]
    extract_GIs()