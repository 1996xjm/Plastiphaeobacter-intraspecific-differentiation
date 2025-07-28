"""
# NCBI染色体质粒命名规则
https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#chr_names
Plasmid names
 - Can contain only digits, dots, underscores, and ASCII characters in plain text in the standard English alphabet.
 - Should start with lower case 'p' UNLESS the plasmid name is not known. In that case use 'unnamed', or "unnamed1" & "unnamed2" for distinct unnamed plasmids.
 - Cannot include the word 'plasmid'
 - Are limited to 20 characters
"""
from Bio import SeqIO
import os

from os import path


def rename(genome_path,  fasta_output_path):

    os.makedirs(fasta_output_path,exist_ok=True)

    for genome in genome_path:

        records = SeqIO.parse(genome, "fasta")

        faa_dict = {}

        for rec in records:
            if "plasmid" in rec.id and "*" not in rec.seq:
                p_id = "_".join(rec.id.split("|")[1].split("_")[:4])
                if p_id in faa_dict:
                    faa_dict[p_id].append(rec)
                else:
                    faa_dict[p_id] = [rec]

        for key, value in faa_dict.items():
            with open(fasta_output_path + f"/{key}.fasta", "w") as f:
                SeqIO.write(value, f, "fasta")


def gbk_chromosome(gbk_path,  fasta_output_path):

    os.makedirs(fasta_output_path,exist_ok=True)

    for gbk in gbk_path:

        records = SeqIO.parse(gbk, "genbank")
        for record in records:
            # print("ID:", record.id)
            # print("Name:", record.name)
            # print("Description:", record.description)
            # print("Number of features:", len(record.features))
            # print("Sequence length:", len(record.seq))
            if "chromosome" in record.id:
                with open(fasta_output_path + f"/{record.id}.gbk", "w") as f:
                    SeqIO.write([record], f, "genbank")



def gbk_plasmid(gbk_path,  fasta_output_path):



    for gbk in gbk_path:

        sample = gbk.split("/")[1]

        out_dir = path.join(fasta_output_path, sample)
        os.makedirs(out_dir, exist_ok=True)

        records = SeqIO.parse(gbk, "genbank")
        for record in records:
            # print("Description:", record.description)
            # print("ID:", record.id)
            if "plasmid" in record.id:
                with open(path.join(out_dir, f"{record.id}.gbk"), "w") as f:
                    SeqIO.write([record], f, "genbank")



if __name__ == '__main__':

    rename(snakemake.input["translated_faa"],  snakemake.output["faa"])
    gbk_chromosome(snakemake.input["annot_gbk"],  snakemake.output["chromosome_gbk"])
    gbk_plasmid(snakemake.input["annot_gbk"],  snakemake.output["plasmid_gbk"])