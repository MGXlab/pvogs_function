#!/usr/bin/env python

import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Subset the scores to a given list of interactions')

    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")
    requiredArgs.add_argument('-s', '--scores-file',
                              dest='scores_file',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The output of calculate_all_scores.py"
                              )
    requiredArgs.add_argument('-i', '--input-ints',
                              dest='input_ints',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A tsv file containing interaction mappings between ncbi and pvogs")
    requiredArgs.add_argument('-o', '--output-file',
                              dest='outfile',
                              required=True,
                              type=lambda p: Path(p).resolve(),
                              help="File path to write the results in")

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()

def parse_pvogs_interactions(interactions_fp):
    interactions = []
    uniques = 0
    all_in = 0
    with open(interactions_fp, 'r') as fp:
        for line in fp:
            all_in += 1
            fields = line.split('\t')
            interaction = tuple(sorted([fields[1].strip(), fields[3].strip()]))
            if interaction[0] == interaction[1]:
                print("Skipping: {} (self)".format(interaction))
                
            elif interaction not in interactions:
                interactions.append(interaction)
                uniques +=1
            else:
                print("Skipping: {} (duplicate)".format(interaction))
    print("Parsed {} / {} total input interactions".format(uniques, all_in))
    return interactions

def subset_all_scores(all_scores_fp, interactions_list, interactions_scores_fp):
    counter = 0 
    with open(all_scores_fp, 'r') as fin, open(interactions_scores_fp, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            fields = line.split('\t')
            interaction = tuple(sorted([fields[0].strip(), fields[1].strip()]))
            if interaction in interactions_list:
                fout.write(line)
                counter += 1
    return counter

def main():
    args = parse_args()
    ints_list = parse_pvogs_interactions(args.input_ints)
    total = subset_all_scores(args.scores_file, ints_list, args.outfile)
    print("{} interaction scores written to file {}".format(total, str(args.outfile)))

if __name__ == '__main__':
    main()

