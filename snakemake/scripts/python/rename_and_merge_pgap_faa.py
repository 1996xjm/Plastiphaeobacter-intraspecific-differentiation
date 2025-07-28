from Bio import SeqIO
from log import Log

def transform_faa(faa_file_list, out_file_faa):

    all_records_list = []

    for faa_file in faa_file_list:
        records = SeqIO.parse(faa_file, "fasta")

        for record in records:
            new_id = "_".join(record.id.split("_prot_")[1].split("_")[:-1])
            record.id = new_id
            # record.description = ""
            all_records_list.append(record)


    with open(out_file_faa, "w") as f:
        SeqIO.write(all_records_list, f, "fasta")

    Log.info(f"total aas count: {len(all_records_list)}")
def transform_fna(fna_file_list, out_file_fna):

    all_records_list = []

    for faa_file in fna_file_list:
        records = SeqIO.parse(faa_file, "fasta")

        for record in records:
            print(record.id)
            new_id = "_".join(record.id.split("_cds_")[1].split("_")[:-1])
            record.id = new_id
            # record.description = ""
            all_records_list.append(record)


    with open(out_file_fna, "w") as f:
        SeqIO.write(all_records_list, f, "fasta")

    Log.info(f"total fnas count: {len(all_records_list)}")


if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        sys.stdout = f
        sys.stderr = f
        transform_faa(snakemake.input["faa"], snakemake.output["faa"])
        if "fna" in snakemake.input:
            transform_fna(snakemake.input["fna"], snakemake.output["fna"])