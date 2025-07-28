
from Bio import SeqIO
import datetime
from log import Log
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def concatenate_seq(OG_file_list, outputFasta):

    record_list = []
    new_record_list = []
    id_list = []
    for i in OG_file_list:
        OG_records_dict = SeqIO.index(i, "fasta")
        record_list.append(OG_records_dict)
        id_list.extend(list(OG_records_dict.keys()))
    id_list = list(set(id_list))

    print(id_list)

    for id in id_list:
        new_seq = Seq("")
        for record in record_list:
            if id in record:
                new_seq += record[id].seq
            else:
                new_seq += Seq("-"*len(record[list(record.keys())[0]].seq))
        new_record = SeqRecord(
            new_seq,
            id=id,
            description=""
        )
        new_record_list.append(new_record)

    SeqIO.write(new_record_list, outputFasta, "fasta")



if __name__ == '__main__':
    concatenate_seq(snakemake.input, snakemake.output[0])