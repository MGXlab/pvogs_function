#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Produce a list of uniprot ids from the metadata.tsv (output of\n"
                                             "summarize_intact.py",
                                formatter_class=argparse.RawTextHelpFormatter)
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--input_tsv',
                    dest='input_tsv',
                    type=str,
                    required=True,
                    help="The metadata.tsv file as produced from summarize_intact.py")
requiredArgs.add_argument('-l', '--output-list',
                    dest='output_list',
                    help="A txt file containing the raw uniprot ids from intact.\n"
                         "This will be used to post a query to the uniprot mapping tool",
                    required=True)
requiredArgs.add_argument('-o', '--output-tsv',
                          dest='output_tsv',
                          required=True,
                          help='A tsv file that contains the result of filtering')

parser._action_groups.append(optionalArgs)

if __name__ == '__main__':
    args = parser.parse_args()
    df = pd.read_csv(args.input_tsv, sep='\t')
    # Select non-self interactions where both interactors come from uniprot
    interactions = df.loc[(df['same_protein'] == 0) &
                          (df['source_A'] == 'uniprotkb') &
                          (df['source_B'] == 'uniprotkb'), ]

    interactions.to_csv(args.output_tsv, sep='\t', index=False)

    print(f'Number of non-self interactions with uniprotkb identifiers: {interactions.shape[0]}')

    # Create a protA_info data-frame with two columns prot_id, source_db
    protA_info = interactions[['prot_A', 'source_A']].reset_index(drop=True)
    protA_info = protA_info.rename(columns={'prot_A': 'prot_id',
                                            'source_A': 'source_db'})
    # Create a protB_info data-frame with two columns prot_id, source_db
    protB_info = interactions[['prot_B', 'source_B']].reset_index(drop=True)
    protB_info = protB_info.rename(columns={'prot_B': 'prot_id',
                                            'source_B': 'source_db'})
    # Concatenate the two data-frames
    prot_info = pd.concat([protA_info, protB_info], ignore_index=True)
    # Drop the duplicates based on uniprot_id
    prot_info.drop_duplicates(subset='prot_id', inplace=True)
    print(f'{prot_info.shape[0]} unique proteins come from uniprot')

    uniprot_ids = prot_info['prot_id'].to_list()

    # Write the identifiers to the output file
    with open(args.output_list, 'w') as fout:
        for uniprot_id in uniprot_ids:
            fout.write(f'{uniprot_id}\n')
