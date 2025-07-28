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


def rename(genome_path, plasmid_dir, chromosome_dir):

    os.makedirs(plasmid_dir, exist_ok=True)
    os.makedirs(chromosome_dir, exist_ok=True)

    for genome in genome_path:

        records = SeqIO.parse(genome, "fasta")

        for rec in records:
            if "plasmid" in rec.id:
                with open(plasmid_dir + f"/{rec.id}.fasta", "w") as f:
                    SeqIO.write([rec], f, "fasta")

            if "chromosome" in rec.id:
                with open(chromosome_dir + f"/{rec.id}.fasta", "w") as f:
                    SeqIO.write([rec], f, "fasta")







if __name__ == '__main__':

    rename(snakemake.input,  snakemake.output["plasmid_dir"], snakemake.output["chromosome_dir"])