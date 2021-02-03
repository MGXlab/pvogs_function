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


if __name__ == "__main__":
    args = parser.parse_args()
    input_fp = Path(args.input_tsv)
    square_matrix_fp = Path(args.output_tsv)
    process_comparem_output(input_fp, square_matrix_fp)

