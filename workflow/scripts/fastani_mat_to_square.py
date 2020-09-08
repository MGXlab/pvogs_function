#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Convert the matrix output of fastANI (v. 1.3) to a square matrix')

    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")
    requiredArgs.add_argument('-i', '--input-matrix',
                              dest='mat_in',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A matrix file from fastANI, produced with its -matrix option")
    requiredArgs.add_argument('-o', '--output-matrix',
                              dest='mat_out',
                              type=lambda p: Path(p).resolve(),
                              required=True,
                              help="File name of the square matrix")
    optionalArgs.add_argument('--process-names',
                              action='store_true',
                              dest='process_names',
                              help='Use the accession as the row/column names and not the full path'
                                   'Assumes that the results look like /path/to/<accession>.fasta'
                              )

    parser._action_groups.append(optionalArgs)
    return parser.parse_args()

def convert_ani_matrix_to_square(ani_output_fp, process_names=True):
    """
    Helper function to convert the output fastANI matrix to a square matrix.
    Makes some calculations easier.
    :param ani_output_fp: Path to fastANI output matrix file
    :return: A tuple of (the matrix itself, the names of the files)
    """
    with open(ani_output_fp, 'r') as fin:
        # first line is the number of genomes compared
        no_of_elements = int(fin.readline().strip())
        # The rest of the lines contain the values
        data_lines = [line.strip() for line in fin]
        # Each line's values start from the second column
        data_fields = [line.split('\t')[1:] for line in data_lines]
        # First element is the genome_file
        genomes_files = [line.split('\t')[0] for line in data_lines]
        if process_names:
            names = []
            for gf in genomes_files:
                # I am expecting something like 'NC_000866.4.fasta'
                fname = Path(gf).name
                name = fname.rstrip('.fasta')
                names.append(name)
        else:
            names = genomes_files 
        # Initialize a square matrix of shape no_of_elements x no_of_elements
        # filled with zeros
        mat = np.zeros((no_of_elements, no_of_elements))

    # A loop to get values to fill the zero-filled matrix
    for i in range(0, no_of_elements):
        for j in range(0, no_of_elements):
            # Diagonal
            if i == j:
                value = 100.
            # Weird but works
            # TO DO add better documentation of what is happening here
            if i < j:
                value = data_fields[j][i]
                if value == 'NA':
                    value = 0.
            if i > j:
                value = data_fields[i][j]
                if value == 'NA':
                    value = 0.
            mat[i, j] = value
    return mat, names

def rewrite_ani_matrix_as_square(fastani_raw, output_fp, **kwargs):
    """
    Helper function for writing the result
    """
    # Get the matrix and names to use
    mat, names = convert_ani_matrix_to_square(fastani_raw, process_names=kwargs['process_names'])
    # Create a data frame
    df = pd.DataFrame(mat, index=names, columns=names)
    # Write the matrix to the file
    df.to_csv(output_fp, header=True, index=True, sep='\t')

def main():
    args = parse_args()
    rewrite_ani_matrix_as_square(args.mat_in, args.mat_out, process_names=args.process_names)

if __name__ == '__main__':
    main()
