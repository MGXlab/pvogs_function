#!/usr/bin/env python

import argparse
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import time
from math import ceil

parser = argparse.ArgumentParser(description='Download a list of ncbi accessions to the output file')
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--input-list',
                          dest='input_list',
                          required=True,
                          help="A txt file containing accessions to get from genbank")
requiredArgs.add_argument('-o', '--output-file',
                    dest='output_file',
                    type=str,
                    required=True,
                    help="The output file to write sequences to")
optionalArgs.add_argument('-e', '--e-mail',
                    help='E-mail address to be used with Bio.Entrez.email. Required by NCBI to notify if something is '
                         'off or if you overload their servers',
                    dest='email',
                    type=str,
                          )
optionalArgs.add_argument('--output-fmt',
                         required=False,
                         dest='output_fmt',
                         default='gb',
                         type=str,
                         help="Store the results in this file format [default='gb' (genbank)]"
                          )

parser._action_groups.append(optionalArgs)

def sequence_info_to_dic(sequence_info_file):
    """
    Collapse protein ids per genome
    Args:
        sequence_info_file (str): A tsv file containing 3 columns uniprot_id, genome_id, protein_id,
                                    with a header.
    Return:
        sequence_info (dict): A dictionary of the form {genome_id: [(protein_id_1, uniprot_id_1),
                                                                    ... ],
                                                                    ... }
    """
    sequence_info = {}
    with open(sequence_info_file, 'r') as f:
        # Skip header
        next(f)
        for line in f:
            fields = [field.strip() for field in line.split('\t')]
            uniprot_id = fields[0]
            genome_id = fields[1]
            protein_id = fields[2]
            if genome_id not in sequence_info:
                sequence_info[genome_id] = [(uniprot_id, protein_id,)]
            else:
                sequence_info[genome_id].append((uniprot_id, protein_id))
    return sequence_info

def txt_file_to_list(genomes_txt):
    """
    Read ids from a one column file to a list

    Args:
        genomes_txt:str: Path to file with the ids.
    """
    with open(genomes_txt, 'r') as fin:
        genomes_list = [line.strip() for line in fin]
    return genomes_list


def download_sequences(genomes_list,
                       genomes_file,
                       email_address='',
                       output_fmt="gb",
                       batch_size = 100
                       ):

    # Required by Bio.Entrez
    if email_address:
        Entrez.email = email_address

    # Some progress tracking
    total_batches = ceil(len(genomes_list) / batch_size)
    batch_no = 0

    with open(genomes_file, 'w') as fout:
        for i in range(0, len(genomes_list), 100):
            batch_no += 1
            batch = genomes_list[i:i + 100]
            print('Downloading batch {}/{}'.format(batch_no, total_batches))
            handle = Entrez.efetch(db="nuccore", id=batch, rettype=output_fmt, retmode="text")
            batch_data = handle.read()
            fout.write(batch_data)
            handle.close()
            # Wait 2 seconds before next batch
            # This is ok for small sets of batches
            time.sleep(2)

def main():
    args = parser.parse_args()

    genomes_list = txt_file_to_list(args.input_list)
    download_sequences(genomes_list, 
                       args.output_file, 
                       args.email, 
                       output_fmt=args.output_fmt)


if __name__ == '__main__':
    main()
