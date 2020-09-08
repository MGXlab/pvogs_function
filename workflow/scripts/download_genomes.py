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
requiredArgs.add_argument('-e', '--e-mail',
                    help='E-mail address to be used with Bio.Entrez.email. Required by NCBI to notify if something is '
                         'off or if you overload their servers',
                    dest='email',
                    type=str,
                    required=True
                    )
optionalArgs.add_argument('--output-fmt',
                         required=False,
                         dest='output_fmt',
                         default='gb',
                         type=str,
                         help="Store the results in this file format [default='gb' (genbank)]")

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
    with open(genomes_txt, 'r') as fin:
        genomes_list = [line.strip() for line in fin]
    return genomes_list


def download_sequences(genomes_list, genomes_file, email_address, output_fmt="gb", batch_size = 100):

    # Required by Bio.Entrez
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
            time.sleep(2)


if __name__ == '__main__':
    args = parser.parse_args()

    genomes_list = txt_file_to_list(args.input_list)
    download_sequences(genomes_list, args.output_file, args.email, output_fmt=args.output_fmt)

    # genomes_fna = args.prefix_out + '_genomes.fna'
    # proteins_faa = args.prefix_out + '_proteins.faa'
    # metadata = args.prefix_out + '_metadata.tsv'
    #
    # with open(genbank_file, 'r') as fin, \
    #         open(genomes_fna, 'w') as fg, \
    #         open(proteins_faa, 'w') as fp, \
    #         open(metadata, 'w') as fout:
    #     # Write the header to the metadata file
    #     fout.write(metadata_header + '\n')
    #
    #     # Loop through the genbank records to extract the information
    #     for record in SeqIO.parse(fin, format="genbank"):
    #         unversioned = record.id.split('.')[0]  # some entries don't come with their version
    #         # Check if the versioned accession is in the dictionary
    #         if record.id in sequence_info:
    #             record_id = record.id
    #         else:  # otherwise use the unversioned string
    #             record_id = unversioned
    #
    #         # Get a list of associated proteins with the record
    #         proteins = [prot_info[1] for prot_info in sequence_info.get(record_id)]
    #
    #         # Create a mapping of protein INSDC accession to uniprot accession
    #         prot_map = {prot_info[1]: prot_info[0] for prot_info in sequence_info.get(record_id)}
    #
    #         # Polyproteins come with one accession
    #         # TO DO
    #         # Get the specific location of each chain
    #         # after it is post-translationally modified
    #         # Maybe  use the FT attribute in the uniprot file?
    #         if len(proteins) > 1 and len(set(proteins)) == 1:
    #             # For now print these out and skip
    #             print("Polyprotein? Nucleotide accession: {}, Protein list: {}".format(record_id, proteins))
    #             pass
    #
    #         else:
    #             # Get the genome sequence anyway
    #             SeqIO.write(record, fg, "fasta")
    #             # The list of protein accessions associated with the genome
    #             # TO DO
    #             # turning them in a set first otherwise there are some duplicates
    #             lproteins = list(set(proteins))
    #             for p in lproteins:
    #                 # The dictionary of features of the genbank record
    #                 for f in record.features:
    #                     # If there is a protein_id field in the qualifiers
    #                     # and is the same as the protein
    #                     if ('protein_id' in f.qualifiers) and (p in f.qualifiers.get('protein_id')):
    #                         # Extract genome length, start, end, strand of the feature
    #                         str_out = '\t'.join(map(str, [p, prot_map[p], record.id, len(record.seq),
    #                                                       f.location.start.position,
    #                                                       f.location.end.position,
    #                                                       f.strand]))
    #                         fout.write(str_out + '\n')
    #                         # Get the translation of the feature
    #                         prot_seq = f.qualifiers.get('translation')[0]
    #                         # Construct a SeqRecord object
    #                         prot_rec = SeqRecord(Seq(prot_seq, IUPAC.protein),
    #                                              id=p,
    #                                              description='')
    #                         # Write the sequence to the proteins file
    #                         SeqIO.write(prot_rec, fp, "fasta")
