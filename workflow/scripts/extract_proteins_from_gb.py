#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


parser = argparse.ArgumentParser(description='Get all protein sequences from a genbank file with'
                                             'annotated genomes.')
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--input-gb',
                          dest='input_gb',
                          required=True,
                          help="A genbank file containing annotated genomes")
requiredArgs.add_argument('-o', '--output-fasta',
                    dest='output_fasta',
                    type=str,
                    required=True,
                    help="The output file to write sequences to")

parser._action_groups.append(optionalArgs)


def extract_all_proteins_from_genbank_file(genbank_fin, fasta_out):
    genomes, proteins = 0, 0
    with open(fasta_out, 'w') as fout:
        for record in SeqIO.parse(genbank_fin, "genbank"):
            genomes += 1
            for f in record.features:
                if f.type == 'CDS':
                    protein = f.qualifiers.get('protein_id')
                    if protein is not None:
                        prot_seq = f.qualifiers.get('translation')[0]
                        prot_rec = SeqRecord(Seq(prot_seq, IUPAC.protein),
                                             id=protein[0],
                                             description='')
                        proteins += SeqIO.write(prot_rec, fout, "fasta")

    return genomes, proteins


if __name__ == '__main__':
    args = parser.parse_args()
    g, p = extract_all_proteins_from_genbank_file(args.input_gb, args.output_fasta)
    print("Extracted {} from {} input genomes".format(p, g))
