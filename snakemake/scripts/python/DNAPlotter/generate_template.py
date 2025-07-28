from Bio import SeqIO
from datetime import datetime
from os import path



def generate():
    records = SeqIO.parse(gbk_file, "genbank")
    for record in records:
        seq_len = len(record.seq)


    if seq_type == "chromosome":

        tmp_str = f"""## DNA Plot :: track template (created: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")})
# line attributes: start=0 end={seq_len} line_size=1 circular=true
# tick marks: major=400000 minor=200000
# GC Graph: height=0.2 window_size=10000 base_step_size=200 track=0.22 minus_colour=255:217:47 plus_colour=166:216:84
# Columns are:
# POS  - track position
# SIZE - track size
# FWD  - show forward strand features
# REV  - show reverse strand features
# NOT  - use NOT
# ANY  - show any features
# KEY  - show features of this key
# QUAL - show features with this qualifier
# VAL  - show features with this qualifier value(s)
# COL  - colour for this track e.g. 255:0:0 (R:G:B) or NULL
# NAME - file entry name or null
# DIR  - root directory for this file
#
#POS	SIZE	FWD	REV	NOT 	ANY	KEY	QUAL	VAL	COL	NAME	DIR
0.95	10.0	true	false	false	false	CDS	null		102:194:165	{path.basename(gbk_file)}	{path.dirname(gbk_file)}
0.9	10.0	false	true	false	false	CDS	null		252:141:98	{path.basename(gbk_file)}	{path.dirname(gbk_file)}
0.85	10.0	true	true	false	false	CDS	null		106:61:154	GIs.gff	{path.abspath(gff_dir)}
0.8	10.0	true	true	false	false	CDS	null		88:187:204	chromosome_GTA.gff	{path.abspath(gff_dir)}
0.8	10.0	true	true	false	false	CDS	null		129:199:132	chromosome_phage_regions_intact.gff	{path.abspath(gff_dir)}
0.8	10.0	true	true	false	false	CDS	null		202:198:0	chromosome_phage_regions_questionable.gff	{path.abspath(gff_dir)}
0.8	10.0	true	true	false	false	CDS	null		255:90:90	chromosome_phage_regions_incomplete.gff	{path.abspath(gff_dir)}
0.75	10.0	true	true	false	false	CDS	null		120:0:1	transposase_genes.gff	{path.abspath(gff_dir)}
0.75	10.0	true	true	false	false	CDS	null		228:26:28	recombinase_integrase_genes.gff	{path.abspath(gff_dir)}
"""
    else:
        tmp_str = f"""## DNA Plot :: track template (created: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")})
# line attributes: start=0 end={seq_len} line_size=1 circular=true
# tick marks: major=10000 minor=5000
# GC Graph: height=0.2 window_size=10000 base_step_size=200 track=0.22 minus_colour=255:217:47 plus_colour=166:216:84
# Columns are:
# POS  - track position
# SIZE - track size
# FWD  - show forward strand features
# REV  - show reverse strand features
# NOT  - use NOT
# ANY  - show any features
# KEY  - show features of this key
# QUAL - show features with this qualifier
# VAL  - show features with this qualifier value(s)
# COL  - colour for this track e.g. 255:0:0 (R:G:B) or NULL
# NAME - file entry name or null
# DIR  - root directory for this file
#
#POS	SIZE	FWD	REV	NOT 	ANY	KEY	QUAL	VAL	COL	NAME	DIR
0.95	10.0	true	false	false	false	CDS	null		102:194:165	{path.basename(gbk_file)}	{path.dirname(gbk_file)}
0.9	10.0	false	true	false	false	CDS	null		252:141:98	{path.basename(gbk_file)}	{path.dirname(gbk_file)}
0.85	10.0	true	true	false	false	CDS	null		120:0:1	transposase_genes.gff	{path.abspath(gff_dir)}
0.85	10.0	true	true	false	false	CDS	null		228:26:28	recombinase_integrase_genes.gff	{path.abspath(gff_dir)}
"""
    with open(orthologous_gene_track_template_file, "r") as f:
        tmp_str = tmp_str + f.read()
    with open(tmp_out, "w") as f:
        f.write(tmp_str)




if __name__ == '__main__':
    gbk_file = snakemake.input["gbk"]
    orthologous_gene_track_template_file = snakemake.input["orthologous_gene_track_template"]
    gff_dir = snakemake.params["gff_dir"]
    tmp_out = snakemake.output[0]
    ref_sample = snakemake.wildcards["sample"]
    seq_type = snakemake.params["seq_type"]


    generate()