from Bio import SeqIO
import pandas as pd
from log import Log
import os




# chromosomal起始点基因：dnaA (GO:0003677)
# plasmid选取起始点基因优先顺序：repA -> parA -> 第一个基因间隔区

def fix(input, output, logPath):
    log_f = open(logPath, "w")
    df_gff = pd.read_table(input[1], skiprows=5)

    df_gff = df_gff.groupby(by="#Sequence Id")

    records = SeqIO.parse(input[0], "fasta")
    records_list = []
    for rec in records:

        topology_description = rec.description.split(" ")[-1]

        if "topology" not in topology_description:
            raise ValueError("Can not find topology description!")


        if "circular" in topology_description:
            Log.info(f"fix {rec.id}.", log_f)
            contig_id = rec.id
            group = df_gff.get_group(contig_id)

            if "location=chromosome" in rec.description:
                # dnaA和dnaN一般相邻出现co-localize
                dnaN_df = group[(group["Gene"] == "dnaN")]
                if len(dnaN_df) > 1:
                    raise ValueError("More than one dnaA gene were found!")
                dnaN_index = dnaN_df.index[0]
                dnaA_df = group[(group["Gene"] == "dnaA")]
                if (dnaN_index-1) in dnaA_df.index:
                    dnaA_row = dnaA_df.loc[dnaN_index - 1, :]
                elif (dnaN_index+1) in dnaA_df.index:
                    dnaA_row = dnaA_df.loc[dnaN_index + 1, :]
                else:
                    raise ValueError("Can not find dnaA gene in chromosome!")



                Log.info(f"Found dnaA: \n"+dnaA_row.to_csv(sep="\t", index=False), log_f)
                start = dnaA_row["Start"]
                stop = dnaA_row["Stop"]
                if dnaA_row["Strand"] == "-":

                    rec.seq = rec.seq[stop:] + rec.seq[:stop]
                    rec.seq = rec.seq.reverse_complement()

                else:
                    rec.seq = rec.seq[start-1:] + rec.seq[:start-1]

            else:
                repA_row = group[(group["Gene"] == "repA")]
                parA_row = group[(group["Gene"] == "parA")]

                if len(repA_row) != 0 or len(parA_row) != 0:
                    if len(repA_row) != 0 :
                        gene_row = repA_row
                    elif len(parA_row) != 0 :
                        gene_row = parA_row


                    Log.info(f"Found repA/parA: \n" + gene_row.to_csv(sep="\t", index=False), log_f)


                    start = gene_row.iloc[0, 2]
                    stop = gene_row.iloc[0, 3]
                    if list(gene_row["Strand"])[0] == "-":

                        rec.seq = rec.seq[stop:] + rec.seq[:stop]
                        rec.seq = rec.seq.reverse_complement()

                    else:
                        rec.seq = rec.seq[start-1:] + rec.seq[:start-1]
                else:
                    Log.warning("Can not find repA or parA gene in plasmid!", log_f)
                    row_index = 0
                    # 寻找第一个间隔区
                    while row_index < group.shape[0]-1:
                        if group.iloc[row_index, 3] < group.iloc[row_index + 1, 2]:
                            end_gene_info = group.iloc[row_index].to_csv(sep='\t')#间隔区前面的基因
                            Log.info(f"found end gene:\n######\n{end_gene_info}######", log_f)
                            end = group.iloc[row_index, 3]
                            rec.seq = rec.seq[end:] + rec.seq[0:end]
                            break
                        row_index += 1

                    if row_index == group.shape[0]-1:
                        Log.warning("There is no end gene!", log_f)



        else:
            Log.warning(f"{rec.id}: This seq is not circular!", log_f)

        records_list.append(rec)
    with open(output, "w") as f:
        SeqIO.write(records_list, f, "fasta")
    log_f.close()


if __name__ == '__main__':
    # input 不支持字典 还是列表
    fix(snakemake.input, snakemake.output[0], snakemake.log[0])