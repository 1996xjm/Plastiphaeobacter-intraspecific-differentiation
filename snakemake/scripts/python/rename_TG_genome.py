"""
# NCBI染色体质粒命名规则
https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#chr_names
Plasmid names
 - Can contain only digits, dots, underscores, and ASCII characters in plain text in the standard English alphabet.
 - Should start with lower case 'p' UNLESS the plasmid name is not known. In that case use 'unnamed', or "unnamed1" & "unnamed2" for distinct unnamed plasmids.
 - Cannot include the word 'plasmid'
 - Are limited to 20 characters
"""
import sys
from Bio import SeqIO
import pandas as pd
from log import Log


def rename(genome_path, assembly_info_path, sample_name, polish_times, fasta_output_path, bakta_replicons_output_path):
    # print(input.genome)
    df = pd.read_table(assembly_info_path)
    records = SeqIO.index(genome_path, "fasta")
    plasmid_index = 1
    contig_index = 1
    records_list = []

    # bakta replicons table

    bakta_replicons_list = []


    for i, row in df.iterrows():

        rec = records.get(f"{row['#seq_name']}{'_pilon' * polish_times}")

        description_list = []



        bakta_seq_type = "contig"

        # assembly_info.txt是有按照contig长短排序的，排第一个的是最长的染色体
        # 这里还没有考虑多条染色体的情况
        if i == 0:
            rec.id = f"{sample_name}_chromosome_1"
            bakta_seq_type = "chromosome"
            description_list.append("[location=chromosome]")
            if row["circ."] != "Y":
                print("##### chromosome is not circular!")
        else:
            description_list.append("[location=plasmid]")
            if row["circ."] == "Y":
                rec.id = f"p{sample_name}_{plasmid_index}"
                bakta_seq_type = "plasmid"
                description_list.append(f"[plasmid-name=p{sample_name}_{plasmid_index}]")
                plasmid_index += 1
            else:
                rec.id = f"{sample_name}_contig_{contig_index}"
                description_list.append(f"[plasmid-name=pc_{sample_name}_{contig_index}]")
                contig_index += 1



        topology = "circular" if row["circ."] == "Y" else "linear"

        description_list.append(f"[flye_seq_name={row['#seq_name']}]")
        description_list.append(f"[length={len(rec.seq)}]")
        description_list.append(f"[coverage={row['cov.']}]")
        description_list.append(f"[topology={topology}]")
        rec.description = " ".join(description_list)
        Log.info(rec.description)

        records_list.append(rec)
        bakta_replicons_list.append(pd.Series({
            "original locus id": rec.id,
            "new locus id": rec.id,
            "type": bakta_seq_type,
            "topology": topology,
            "name": "-"
        }))
    bakta_replicons_df = pd.concat(bakta_replicons_list, axis=1).T.set_index("original locus id", drop=True)

    with open(fasta_output_path, "w") as f:
        SeqIO.write(records_list, f, "fasta")

    bakta_replicons_df.to_csv(bakta_replicons_output_path, sep="\t")


if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        sys.stdout = f
        sys.stderr = f
        rename(snakemake.input[0], snakemake.params.assembly_info, snakemake.params.sample, snakemake.params.polish_times, snakemake.output["fasta"], snakemake.output["bakta_replicons_table"])