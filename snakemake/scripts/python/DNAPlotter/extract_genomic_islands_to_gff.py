from os import path
import pandas as pd
import re
def extract_GIs():
    backbone_col_name_list = []
    with open(xmfa_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            line = line.strip()

            if re.match(r'^#Sequence\d+File', line):
                sample_name = "_".join(path.basename(line.split("\t")[1]).split("_")[:2])
                backbone_col_name_list.append(f"{sample_name}_leftend")
                backbone_col_name_list.append(f"{sample_name}_rightend")
    backbone_df = pd.read_table(backbone_file)
    backbone_df.columns = backbone_col_name_list
    backbone_filtered_df = backbone_df.drop(columns=[f"{ref_sample}_leftend", f"{ref_sample}_rightend"])
    GIs_df = backbone_df[(backbone_filtered_df == 0).all(axis=1)][
        [f"{ref_sample}_leftend", f"{ref_sample}_rightend"]]
    GI_index = 1
    GI_gff_f = open(gff_out_path, "w")
    for i,row in GIs_df.iterrows():
        leftend = row[f"{ref_sample}_leftend"]
        rightend = row[f"{ref_sample}_rightend"]
        print("\t".join([f"{ref_sample}_chromosome_1", "Local", "CDS", f"{leftend}", f"{rightend}", ".", "+", ".", f"Name={ref_sample}_GI_{GI_index}"]), file=GI_gff_f)
        GI_index += 1

    GI_gff_f.close()


if __name__ == '__main__':

    ref_sample = snakemake.wildcards["sample"]

    xmfa_file = snakemake.input["xmfa_file"]
    backbone_file = snakemake.input["backbone_file"]
    gff_out_path = snakemake.output[0]
    extract_GIs()