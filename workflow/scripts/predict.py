#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd

import pickle

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split


FEATURES = ['jaccard_score',
            'same_score',
            'inwards_score',
            'outwards_score',
            'avg_distance',
            'mean_ani',
            'mean_aai']

# Random State - for reproducibility
RS = 1

def parse_args():
    parser = argparse.ArgumentParser(description='Predict interactions on target')

    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")
    requiredArgs.add_argument('-m', '--model',
                              dest='model_fp',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A pickle file containing the RF parameters"
                              )
    requiredArgs.add_argument('-t', '--target_tsv',
                              dest='target_tsv',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A tsv file containing feature calculations for the interactions")
    requiredArgs.add_argument('-p', '--positive-set',
                              dest='positives_tsv',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="Positive set")
    requiredArgs.add_argument('-n', '--negative-set',
                              dest='negatives_tsv',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="Negative set with labels")
    requiredArgs.add_argument('-o', '--output-file',
                              dest='outfile',
                              required=True,
                              type=lambda p: Path(p).resolve(),
                              help="File path to write the results in")
    optionalArgs.add_argument('-j', '--jobs',
                              dest='n_jobs',
                              type=int,
                              required=False,
                              default=1,
                              help="Number of parallel jobs to run for classification (default : 1)"
                              )
    parser._action_groups.append(optionalArgs)

    return parser.parse_args()


def read_scores_table(scores_fp, label=None):
    """
    Read a table from a tsv file provided as a path.
    Return a dataframe.

    If label is set, append a column named `label` with
    the provided value
    """
    scores_df = pd.read_csv(scores_fp, sep='\t')
    scores_df['interaction'] = scores_df['pvog1'] + '-' + scores_df['pvog2']
    # Unpack features list because I want interaction prepended.
    scores_df = scores_df[['interaction', *FEATURES]]
    if label is not None:
        scores_df['label'] = label
    return scores_df


def concat_data_frames(pos_df,
                       neg_df,
                       subsample=False,
                       clean=True,
                       is_scaled=True,
                      ):
    """
    Concatenate two dataframes

    subsample:bool 
        Subsample the `neg_df` to a number of
        observations equal to the number of pos_df 
        (balance datasets)

    clean:bool
        Remove observations with a value of `avg_distance` == 100000 or 1.

    is_scaled:bool
        The features have been scaled to a range 0-1

    Return:
    concat_df: pd.DataFrame
        The concatenated data frame
    """

    if clean is True:
        pos_df = remove_ambiguous(pos_df)
        neg_df = remove_ambiguous(neg_df)

    n_positives = pos_df.shape[0]
    n_negatives = neg_df.shape[0]

    # Remove possible duplicate interactions from the negatives
    # This might happen because of the random selection when creating the set
    # Why I also select more negatives to begin with
    neg_df = neg_df.loc[~neg_df['interaction'].isin(pos_df['interaction'])]

    if (n_positives != n_negatives) and (subsample is True):
        neg_df = neg_df.sample(n=n_positives, random_state=1)
    concat_df = pd.concat([pos_df, neg_df]).reset_index(drop=True)

    assert concat_df[concat_df.duplicated(subset=['interaction'])].empty == True, concat_df.loc[concat_df.duplicated(subset=['interaction'], keep=False)]

    return concat_df


def scale_df(input_df):
    """
    Scale all feature values in the data frame to [0-1].
    """
    maxes = input_df[FEATURES].max(axis=0)
    scaled_data = input_df[FEATURES].divide(maxes)
    if 'label' in input_df.columns:
        scaled_df = pd.concat([input_df['interaction'], scaled_data, input_df['label']], axis=1)
    else:
        scaled_df = pd.concat([input_df['interaction'], scaled_data], axis=1)
    return scaled_df


def remove_ambiguous(input_df):
    """
    Select observations in the `input_df` that have feature values
    """
    df_clean = input_df[(input_df.jaccard_score != 0)] # This is true if they don't co-occur
    return df_clean



if __name__ == '__main__':
    args=parse_args()

    # Create the training set dataframe
    pos_df = read_scores_table(args.positives_tsv, label=1)
    print("Raw positive set interactions : {}".format(pos_df.shape[0]))

    neg_df = read_scores_table(args.negatives_tsv, label=0)
    print("Raw negative set interactions : {}".format(neg_df.shape[0]))

    scaled_pos_df = scale_df(pos_df)
    scaled_neg_df = scale_df(neg_df)
    training_df = concat_data_frames(scaled_pos_df,
                                     scaled_neg_df,
                                     subsample=True,
                                     clean=True,
                                     is_scaled=True)
    print("Processed training set interactions: {}".format(training_df.shape[0]))

    # Read in the target dataset
    target_df = read_scores_table(args.target_tsv)
    input_targets = target_df.shape[0]
    print("Target set interactions: {} ".format(input_targets))

    # Clean the target from distances
    target_df = remove_ambiguous(target_df)
    # Scale it
    target_df = scale_df(target_df)

    # Remove posnegs
    target_df = target_df.loc[~target_df['interaction'].isin(training_df['interaction'])]
    final_targets = target_df.shape[0]
    print("Removed {} interactions from target ( {} remaining)".format((input_targets - final_targets), final_targets))

    # Read in the model
    with open(args.model_fp, 'rb') as fin:
        RF = pickle.load(fin)

    if args.n_jobs:
        RF.n_jobs = args.n_jobs

    print("Classifier: {}".format(RF))


    Xt = training_df[FEATURES]
    yt = training_df['label']
    # Keep these to append 
    interactions = training_df['interaction']
    # Get the feature values only, for predictions on the target
    X = target_df[FEATURES]

    # Split the training set to train/holdout (0.7/0.3)
    # Holdout is thrown away here...
    # This is required for consistency with the whole model selection process
    X_train, X_holdout, y_train, y_holdout = train_test_split(Xt,
                                                              yt,
                                                              test_size=0.3,
                                                              random_state=RS)
    # write the actual training set used to a file
    # So this can be appended to the complete final table
    # Make a copy
    X_training_out = X_train.copy()
    # Append the label
    X_training_out['label'] = y_train
    # Append the interaction information
    X_training_out['interaction'] = training_df['interaction']
    # Give the training set a proba of 1
    X_training_out['proba'] = 1.
    # Rearrange columns
    X_training_out = X_training_out[['interaction', *FEATURES, 'label', 'proba']]
    
    # Construct the name of the output file
    dir_name = args.outfile.parent
    training_set_fp = dir_name / Path("final_training_set.tsv")
    X_training_out.to_csv(training_set_fp, sep = '\t', index=False)




    print("Fitting...")
    RF.fit(X_train, y_train)

    print("Predicting...")
    X_pred = RF.predict(X)
    X_pred_proba = RF.predict_proba(X)

    target_df['label'] = X_pred
    target_df['proba'] = X_pred_proba[:,1]

    target_df.to_csv(args.outfile, sep="\t", index=False)
    print("Finished! Results are written in {}".format(str(args.outfile)))

