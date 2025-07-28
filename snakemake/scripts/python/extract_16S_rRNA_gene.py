from os import path
from Bio import SeqIO
from log import Log
def extract(rRNA_fasta_list, output_fasta, logPath):
    log_f = open(logPath, "w")
    all_rec = []
    for rRNA_fasta in rRNA_fasta_list:
        sampleName = ".".join(path.basename(rRNA_fasta).split(".")[:-1])
        Log.info(f"#### {sampleName} ####", log_f)
        records = SeqIO.parse(rRNA_fasta, "fasta")

        rRNA_16S_record = []

        flag = False

        for r in records:
            if r.id.startswith("16S_rRNA"):
                id = r.id

                r.id = sampleName
                rRNA_16S_record.append(r)
                flag = True


        if flag:
            # 取第一条16S
            target_index = 0

            if sampleName == "PT24_64":
                target_index = 1


            all_rec.append(rRNA_16S_record[target_index])
            Log.info("Extract successfully!", log_f)
        else:
            Log.warning("There is no 16S!", log_f)

    with open(output_fasta, "w") as f:
        SeqIO.write(all_rec, f, "fasta")

    log_f.close()





if __name__ == '__main__':
    extract(snakemake.input, snakemake.output[0], snakemake.log[0])