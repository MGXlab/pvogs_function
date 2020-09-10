#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Join the predictions.tsv and annotations.tsv "
                                        "in a final table where all information is available "
                                    )
    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument('-i', '--input-file',
                            dest='input_fp',
                            type=str,
                            required=True,
                            help="The predictions.tsv as produced by predict.py")
    requiredArgs.add_argument('-a', '--annotations-file',
                        dest='ann_fp',
                        type=str,
                        required=True,
                        help="The annotations.tsv as produced by process_annotations.py")
    requiredArgs.add_argument('-o', '--output-file',
                              dest='output_fp',
                              required=True,
                              type=str,
                              help="A file to write the final product of joining the two inputs. "
                              "This contains raw and processed annotation per pVOG pair, along with "
                              "the label (0 for no interaction, 1 for predicted interaction) "
                              "and the probability of the prediction."
                              )

    parser._action_groups.append(optionalArgs)
    return parser.parse_args()


FEATURES = ['jaccard_score',
            'avg_distance',
            'same_score',
            'inwards_score',
            'outwards_score',
            'mean_ani',
            'mean_aai',
            ]

def main():
    args = parse_args()
    # Read in the annotations in dataframe and set the index to
    # the pvog column
    ann_df = pd.read_csv(args.ann_fp, sep='\t').set_index('pvog')

    # Read the predictions into a dataframe
    predictions_df = pd.read_csv(args.input_fp, sep="\t")

    # Sort by probability value
    predictions_df = predictions_df.sort_values(by=['proba'], ascending=False)

    # Create a list for each pvog in the interaction pair
    # so they can be used for individual lookups
    pvog1 = [i.split('-')[0] for i in predictions_df['interaction'].values]
    pvog2 = [i.split('-')[1] for i in predictions_df['interaction'].values]

    # Append these lists in the original 
    predictions_df['pvog1'] = pvog1
    predictions_df['pvog2'] = pvog2

    # Get the annotation values from the annotations dataframe
    pvog1_raw = [ann_df.loc[pvog, 'annotation_raw']
                 for pvog in predictions_df['pvog1'].values]

    pvog1_processed = [ann_df.loc[pvog, 'annotation_processed']
                       for pvog in predictions_df['pvog1'].values]
    pvog2_raw = [ann_df.loc[pvog, 'annotation_raw']
                 for pvog in predictions_df['pvog2'].values]
    pvog2_processed = [ann_df.loc[pvog, 'annotation_processed']
                       for pvog in predictions_df['pvog2'].values]

    # Append these to the predictions_df
    predictions_df['pvog1_annotation_raw'] = pvog1_raw
    predictions_df['pvog1_annotation_processed'] = pvog1_processed
    predictions_df['pvog2_annotation_raw'] = pvog2_raw
    predictions_df['pvog2_annotation_processed'] = pvog2_processed

    # Re-arrange the columns just for fun
    predictions_df = predictions_df[
        ['interaction',
         'label',
         'proba',
         'pvog1_annotation_raw',
         'pvog1_annotation_processed',
         'pvog2_annotation_raw',
         'pvog2_annotation_processed',
         *FEATURES,
          ]
    ]

    predictions_df.to_csv(args.output_fp, sep="\t", index=False)

if __name__ == '__main__':
    main()

