import sys
from os import path
sys.path.append(path.join(path.dirname(__file__), ".."))
from log import Log
def gen_mcl_input(blastall_results, mcl_input_file_path, minbit = 0.5, min_percent_identity = 0):

    Log.info(f"minbit: {minbit}")


    all_ids = set([])

    # mapping for the fields in the blast output
    mapping = [str, str, float, float, int, int, int, int, int, float, float]

    # here we perform an initial pass on the blast results to fill the dict that will hold
    # the bit score for each gene when it was blasted against itself. this dictionary
    # will then be used to calculate the 'minbit' value between two genes, which I learned
    # from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8). ITEP defines minbit as
    # 'bit score between target and query / min(selfbit for query, selbit for target)'. This
    # heuristic approach provides a mean to set a cutoff to eliminate weak matches between
    # two genes. minbit value reaches to 1 for hits between two genes that are almost identical.
    self_bit_scores = {}
    line_no = 1
    for line in open(blastall_results):
        fields = line.strip().split('\t')

        try:
            query_id, subject_id, perc_id, q_covs, q_length, s_length, aln_length, mismatches, gaps, e_val, bit_score = \
                [mapping[i](fields[i]) for i in range(0, len(mapping))]
        except Exception as e:
            raise ValueError("dfdsf")
        line_no += 1
        all_ids.add(query_id)
        all_ids.add(subject_id)

        if query_id == subject_id:
            self_bit_scores[query_id] = bit_score



    Log.info(f"all_gene_ids: {len(all_ids)}")

    # CONTINUE AS IF NOTHING HAPPENED


    mcl_input = open(mcl_input_file_path, 'w')

    line_no = 1
    num_edges_stored = 0
    for line in open(blastall_results):
        fields = line.strip().split('\t')

        query_id, subject_id, perc_id, q_covs, q_length, s_length, aln_length, mismatches, gaps, e_val, bit_score = \
            [mapping[i](fields[i]) for i in range(0, len(mapping))]

        line_no += 1


        #
        # FILTERING BASED ON PERCENT IDENTITY
        #

        if perc_id < min_percent_identity:
            continue

        #
        # FILTERING BASED ON MINBIT
        #
        new_minbit = bit_score / min(self_bit_scores[query_id], self_bit_scores[subject_id])
        if new_minbit < minbit:
            continue

        mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
        num_edges_stored += 1

    Log.info(f"num_edges_stored: {num_edges_stored}")



    mcl_input.close()


    return mcl_input_file_path


if __name__ == '__main__':
    with open(snakemake.log[0], 'w') as f:
        sys.stdout = f
        sys.stderr = f
        gen_mcl_input(snakemake.input[0], snakemake.output[0], snakemake.params.minbit)