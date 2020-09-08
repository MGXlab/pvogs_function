#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Rewrite comparem output as square-matrix and calculate "
                                             "mean aai and nonzero mean aai for the records.")
parser.add_argument("--input-tsv", "-i",
                    required=True,
                    dest="input_tsv",
                    help="The aai_summary.tsv file produced with comparem's aai_wf")
parser.add_argument("--output-tsv", "-o",
                    dest='output_tsv',
                    required=True,
                    help="Output tsv file to write square matrix to")

def process_comparem_output(aai_out, matrix_out):
    with open(aai_out, 'r') as fin:
        # Skip header
        next(fin)
        comparisons = {}
        for line in fin:
            fields = [f.strip() for f in line.split()]

            genomeA, genomeB, mean_aai = fields[0], fields[2], float(fields[5])

            if genomeA not in comparisons:
                comparisons[genomeA] = {genomeA: 100.}
            if genomeB not in comparisons:
                comparisons[genomeB] = {genomeB: 100.}

            comparisons[genomeA][genomeB] = mean_aai
            comparisons[genomeB][genomeA] = mean_aai

    df = pd.DataFrame.from_dict(comparisons)
    # Write the square matrix
    df.to_csv(matrix_out, sep='\t', header=True)
    return True

#    # Convert df to numpy array
#    mat = df.to_numpy()
#    # Calculate mean aai for all comparisons
#    aai = mat.mean()
#    # Calculate non zero mean
#    nonzeros = np.count_nonzero(mat)
#    nonzero_mean_aai = mat.sum() / nonzeros
#
#    return aai, nonzero_mean_aai


if __name__ == "__main__":
    args = parser.parse_args()
    input_fp = Path(args.input_tsv)
    square_matrix_fp = Path(args.output_tsv)
    process_comparem_output(input_fp, square_matrix_fp)

#    if args.write_result:
#        # I assume the aai_summary.tsv is in /path/to/INTERACTION/aai_out/aai/aai_summary.tsv
#        interaction = input_dir.parent.parent.name
#        result_file = Path.joinpath(input_dir, 'aai_means.tsv')
#        with open(result_file, 'w') as fout:
#            fout.write('interaction\tmean_aai\tnonzero_mean_aai\n')
#            fout.write(f'{interaction}\t{aai}\t{nonzero_mean_aai}\n')
#    else:
#        print(f'Mean aai : {aai}')
#        print(f'Non-zero mean aai: {nonzero_mean_aai}')
