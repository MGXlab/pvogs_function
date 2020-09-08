#!/usr/bin/env python

import argparse
from itertools import combinations
import random
from Bio import SeqIO
import pathlib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


parser = argparse.ArgumentParser(description='Create a 2-column tsv file with --number-of combinations of proteins'
                                             'from the same genome, for all genomes that have more than 1 '
                                             'protein. '
                                             'If a 2-column tsv files is provided with --exclude, the set of these '
                                             'interactions will be excluded.')
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--input-gb',
                          dest='input_gb',
                          required=True,
                          help="A genbank file containing annotated genomes")
requiredArgs.add_argument('-o', '--fasta-out',
                          dest='fasta_out',
                          type=str,
                          required=True,
                          help="The output faa file to write to")
requiredArgs.add_argument('-x', '--interactions-out',
                          dest='interactions_out',
                          type=str,
                          required=True,
                          help="The output tsv file to write the interactions in (x for A *x* B)")

requiredArgs.add_argument('-a', '--all-proteins-faa',
                          dest='all_proteins_fasta',
                          help='A fasta file that contains all proteins as extracted from the genbank file',
                          required=True)
optionalArgs.add_argument('--number-of-pairs',
                          required=False,
                          help="select this amount of pairs to output",
                          dest='sample_no',
                          type=int)
optionalArgs.add_argument('--exclude',
                          required=False,
                          type=str,
                          help="A 2-column tsv file that contains interactions to be excluded")
optionalArgs.add_argument('--sample-size',
                          required=False,
                          dest='sample_size',
                          type=int,
                          help="If an --exclude file is present, set the number of interactions"
                               "to sample to N * sample-size, where N=number of interactions in the exclude file",
                          default=1
                          )
optionalArgs.add_argument('--random-seed',
                          required=False,
                          dest='random_seed',
                          type=int,
                          default=46,
                          help="Set the random seed for reproducibility [default = 46]"
                         )

parser._action_groups.append(optionalArgs)

def create_all_possible_combos(genbank_fin):
    """
    Create a set of all possible protein interactions, per genome,
    from the input genbank_file
    """
    # all_pairs_dict = {}
    all_pairs = []
    for record in SeqIO.parse(genbank_fin, format='genbank'):
        proteins = []
        for f in record.features:
            if f.type == 'CDS':
                protein = f.qualifiers.get('protein_id')

                if protein:
                    protein_id = protein[0]
                    proteins.append(protein_id)

        if len(proteins) >= 2:
            protein_combos = list(combinations(sorted(proteins), 2))

            # all_pairs_dict[genome_acc] = protein_combos
            for pair in protein_combos:
                all_pairs.append(pair)
        elif len(proteins) == 1:
            pass

    # Sort the tuples by protein id, to make comparisons easier
    all_pairs_sorted = sorted([tuple(sorted(i)) for i in all_pairs])

    return all_pairs_sorted


def create_exclude_set(exclude_file):
    """
    Create a list of interactions to be excluded from analysis,
    e.g. if they are part of the positive set.
    """
    excludes = []
    with open(exclude_file, 'r') as fin:
        for line in fin:
            fields = line.split('\t')
            p1, p2 = fields[0].strip(), fields[1].strip()
            interaction = (p1,p2)
            excludes.append(sorted(interaction))
    exclude_pairs = sorted(set([tuple(sorted(i)) for i in excludes]))

    return exclude_pairs


if __name__ == '__main__':
    args = parser.parse_args()

    interactions_pool = create_all_possible_combos(args.input_gb)
    #print("Total raw pool: {}".format(len(interactions_pool)))

    if args.exclude:
        exclude_set = create_exclude_set(args.exclude)
        print('{} were input for exclusion'.format(len(exclude_set)))
        interactions_pool = sorted(set(interactions_pool) - set(exclude_set))
        print("Total pool after exclusion: {}".format(len(interactions_pool)))
        if not args.sample_size:
            print("The number of interactions to be selected will be set to {}".format(len(exclude_set)))
            print("You can use the --sample-size argument to sample more")
            sample_size = args.sample_size
        else:
            sample_size = args.sample_size * len(exclude_set)
    elif not args.exclude and not args.sample_no:
        # Just testing this
        parser.error("Either provide an input file with interactions to exclude and set the sample-size option"
                     "or set the number of interactions you want to sample with the --number-of-pairs option")
    elif args.sample_no:
        sample_size = args.sample_no
    else:
        parser.error("This is unexpected")

    with open(args.all_proteins_fasta, 'r') as fin:
        seq_data = SeqIO.to_dict(SeqIO.parse(fin, format="fasta"))
    
    random.seed(args.random_seed)

    final_interactions = random.sample(interactions_pool, sample_size)
    final_interactions = sorted([tuple(sorted(i)) for i in final_interactions])
    #print(final_interactions[:5])
    #print("Produced {} interactions".format(len(final_interactions)))

    with open(args.interactions_out, 'w') as fout:
        for interaction in final_interactions:
            fout.write('{}\t{}\n'.format(interaction[0], interaction[1]))
    #print("Interaction data are stored in {}".format(args.interactions_out))

    accessions_list = []
    for interaction in final_interactions:
        accessions_list.append(interaction[0])
        accessions_list.append(interaction[1])
    accessions_list = sorted(set(accessions_list))

    writer = 0
    with open(args.fasta_out, 'w') as fout:
        for accession in accessions_list:
            writer += SeqIO.write(seq_data.get(accession), fout, format="fasta")
    #print("Wrote {} unique sequences in file {}".format(writer, args.fasta_out))
