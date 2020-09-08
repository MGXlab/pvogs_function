#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path
from collections import Counter
import operator

def parse_args():
    parser = argparse.ArgumentParser(description="Parse pvogs annotations from the VOGProteinTable.txt "
                                        "to a table that contains most frequently occurring raw and "
                                        "processed annotation")
    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument('-i', '--input-file',
                            dest='input_fp',
                            type=str,
                            required=True,
                            help="The VOGProteinTable.txt from the database")
    requiredArgs.add_argument('-o', '--output-file',
                        dest='output_fp',
                        type=str,
                        required=True,
                        help="The output table location")

    parser._action_groups.append(optionalArgs)
    return parser.parse_args()

def parse_annotations(vog_table_fp):
    """
    Parse the pvogs annotations from the VOGProteinTable.txt
    This tries to normalize the annotations into a standard
    set of descriptions that is used for counting.
    
    Anything that is annotated as hypothetical is translated
    to 'unknown'.
    
    'protein' and 'putative' are also stripped.
    
    Returns
    annotations_dic:dict
        A dictionary that holds the processed terms
        of the form
            {pvog_id : 
                [annotation_1, annotation_2, ...],
            ...}
            
    annotations_original:dict
        A dictionary that holds the original terms
        of the form
            {pvog_id: 
                [annotation_1, annotation_2, ...],
            ...}
    
    """

    annotations_dic = {}
    annotations_original = {}
    with open(vog_table_fp, 'r') as fin:
        for line in fin:
                fields = line.split('|')
                pvog = fields[0].split(':')[0]
                if len(fields) > 2:
                    annotation = fields[2].split(':')[0]
                    processed_annotation = annotation
                    if 'hypothetical' in processed_annotation:#or annotation == '-':
                        processed_annotation=processed_annotation.replace('hypothetical protein', '').strip()
                        if processed_annotation == '':
                            processed_annotation = 'unknown'
                        else:
                            processed_annotation = processed_annotation.strip()

                    if 'protein' in annotation:
                        processed_annotation = processed_annotation.replace('protein', '')
                    if 'putative' in annotation:
                        processed_annotation = processed_annotation.replace('putative', '')
                else:
                    annotation = 'unknown'
                    processed_annotation = 'unknown'

                if pvog not in annotations_dic:
                    annotations_dic[pvog] = [processed_annotation.strip()]
                else:
                    annotations_dic[pvog].append(processed_annotation.strip())

                if pvog not in annotations_original:
                    annotations_original[pvog] = [annotation.strip()]
                else:
                    annotations_original[pvog].append(annotation.strip())

    return annotations_dic, annotations_original

def get_max_count_annotation(annotations_dic):
    """
    Get the annotation with the max count from a given `annoations_dic`
    
    annotations_dic:dict
        A dictionary of the form 
             {pvog_id : 
                [annotation_1, annotation_2, ...],
            ...}
    Return:
        unique_only:dict
            A dictionary with a pvog as key and
            the annotation with the highest occurring frequency
    """
    unique_only = {}
    for pvog in annotations_dic:
        c = Counter(annotations_dic.get(pvog))
        most_common = max(c.items(), key=operator.itemgetter(1))[0]
        unique_only[pvog] = most_common
    return unique_only

def create_annotations_df(vog_table_fp):
    """
    A wrapper function that parses the annotations file
    into a table

    vog_table_fp:Path-like object
        The filepath of VOGProteinTable.txt

    Return
    df_processed:pd.DataFrame
        A dataframe that holds processed and raw annotations
        for each pvog        
    """
    parsed_annotations = parse_annotations(vog_table_fp)
    processed_annotations = parsed_annotations[0]
    raw_annotations = parsed_annotations[1]
    max_data_processed = get_max_count_annotation(processed_annotations)
    max_data_raw = get_max_count_annotation(raw_annotations)
    df_processed = pd.DataFrame.from_dict(max_data_processed, 
                                          orient='index', 
                                          columns=['annotation_processed'])
    df_raw = pd.DataFrame.from_dict(max_data_raw, 
                                    orient='index', 
                                    columns=['annotation_raw'])
    df_processed = df_processed.reset_index().rename(columns = {'index': 'pvog'})
    df_raw = df_raw.reset_index().rename(columns={'index': 'pvog'})
    df_processed['annotation_raw'] = df_raw['annotation_raw']
    df_processed = df_processed.set_index('pvog')

    return df_processed

def main():
    args = parse_args()
    annotations_df = create_annotations_df(Path(args.input_fp))
    annotations_df.to_csv(Path(args.output_fp), sep='\t')


if __name__ == '__main__':
    main()



