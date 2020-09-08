#!/usr/bin/env python

import argparse
from pathlib import Path
from Bio import SeqIO, SearchIO
import pandas as pd
import numpy as np
from itertools import combinations

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate scores for a set of pVOG interactions, '
                                                  'provided as a tsv file')

    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")
    requiredArgs.add_argument('-p', '--profiles-file',
                              dest='profiles_file',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The all.hmm file from pVOGs database"
                              )
    requiredArgs.add_argument('-g', '--genomes',
                              dest='genomes_fasta',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A fasta file with the genomes")
    requiredArgs.add_argument('-hmm', '--input-hmm',
                              dest='hmmer_in',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The regular output file from hmmsearch all pvogs against"
                                    "the translated genomes"
                             )
    requiredArgs.add_argument('-ani_f', '--ani-matrix',
                              dest='ani_matrix',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The square matrix resulting from fastANI with all genomes"
                             )

    requiredArgs.add_argument('-aai_f', '--aai-matrix',
                              dest='aai_matrix',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The aai square matrix from compareM on all genomes"
                              )
    requiredArgs.add_argument('-o', '--output-file',
                              dest='outfile',
                              required=True,
                              type=lambda p: Path(p).resolve(),
                              help="File path to write the results in")

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()

def get_seq_sizes(seq_fasta):
    """
    Create a dict that holds length sizes for all records
    in a fasta file.
    """
    seq_sizes = {}
    with open(seq_fasta, 'r') as fin:
        for record in SeqIO.parse(fin, "fasta"):
            seq_sizes[record.id] = len(record.seq)
    return seq_sizes

def get_maximum_index(values_list):
    """
    Get the index of the maximum value in a list
    
    param: list: values_list A list of values
    return: int: The index of the maximum value
    """
    if len(values_list) == 1:
        max_index = 0
    else:
        max_index = values_list.index(max(values_list))
    return max_index

def translate_to_genomic_coords(start, end, frame, genome_size):
    """
    Translate the coordinates of the protein from transeq to genomic
    coordinates.
    
    Strand is used here as orientation and not really [non-]conding.
    If the frame is 1,2 or 3 (-->) I call this (+) strand.
    Else, it is the (-) strand
    
    param: int: start The starting coordinate on the protein
    param: int: end The ending coordinate on the protein
    param: int: frame The frame on which it is found (1-6)
    param: int: genome_size The size of the genomes
    
    return: tuple: (genomic start, genomic end, strand)
    """
    nucleic_start = start * 3
    nucleic_end = end * 3
    if frame == 1:
        genomic_start = nucleic_start - 2
        genomic_end = nucleic_end - 2
    if frame == 2:
        genomic_start = nucleic_start - 1
        genomic_end = nucleic_end - 1
    if frame == 3:
        genomic_start = nucleic_start
        genomic_end = nucleic_end
    if frame == 4:
        genomic_start = genome_size - (nucleic_start - 2)
        genomic_end = genome_size -  (nucleic_end - 2)
    if frame == 5:
        genomic_start = genome_size - (nucleic_start - 1)
        genomic_end = genome_size - (nucleic_end -1)
    if frame == 6:
        genomic_start = genome_size - nucleic_start
        genomic_end = genome_size - nucleic_end
        
    if frame in [1,2,3]:
        strand = '+'
    elif frame in [4,5,6]:
        strand = '-'
    else:
        raise ValueError("frame should be one of 1,2,3,4,5,6")
        
    return genomic_start, genomic_end, strand

def collect_hmmsearch_info(hmmer_in, genome_sizes):
    """ 
    Parse the information in a dictionary of the form
    { pvog_id: 
        { genome_id : 
            { frame : [(genomic_start, genomic_end, accuracy, pvog_coverage), ...],
            ... , }
        }
    }
    """
    
    hmm_results = {}
    with open(hmmer_in, 'r') as fin:
        for record in SearchIO.parse(fin, "hmmer3-text"):
            # Initialize an empty dictionary for the pvog
            hmm_results[record.id] = {}
            for hit in record.hits:
                if hit.is_included:
                    # From the transeq output, accessions are suffixed
                    # with _1, _2, _3, ..., depending on the strand
                    genome_id = '_'.join(hit.id.split('_')[0:2])
                    frame = int(hit.id.split('_')[2])
                    hsps_included = [hsp.is_included for hsp in hit.hsps]
                    # For multiple hsps that pass the threshold
                    if any(hsps_included):
                        # Get their score
                        scores = [hsp.bitscore for hsp in hit.hsps]
                        # and select the best one
                        max_i = get_maximum_index(scores)
                        best_hsp = hit.hsps[max_i]
                        
                        # Translate back to genomic coordinates
                        genomic_coords = translate_to_genomic_coords(best_hsp.env_start, 
                                                                     best_hsp.env_end, 
                                                                     frame, 
                                                                     genome_sizes.get(genome_id))
                    
                        span = (best_hsp.env_end - best_hsp.env_start) / record.seq_len                  
                    
                        if genome_id not in hmm_results[record.id]:
                            hmm_results[record.id][genome_id] = (frame, genomic_coords[0],
                                                                 genomic_coords[1], 
                                                                 best_hsp.acc_avg,
                                                                 span, 
                                                                 genomic_coords[2])
    return hmm_results


def get_mean_from_df(subset_list, df):
    df_subset = df.loc[subset_list, subset_list]
    return df_subset.values.mean()


def get_pvogs_ids(profiles_file):
    all_pvogs = []
    with open(profiles_file, 'r') as fin:
        for line in fin:
            if line.startswith('NAME'):
                pvog = line.split()[1]
                all_pvogs.append(pvog)
    return all_pvogs

def get_shortest_distance(startA, endA, startB, endB):
    start_to_start = abs(startA - startB)
    start_to_end = abs(startA - endB)
    end_to_start = abs(endA - startB)
    end_to_end = abs(endA - endB)
    all_distances = [start_to_start, start_to_end, end_to_start, end_to_end]
    return min(all_distances)

def calculate_scores(interaction_tuple, hmm_results, ani_df, aai_df):
    same_score = 0
    inwards_score = 0
    outwards_score = 0
    avg_distance = 1000000 # 1000000 for pairs that are never on the same genome
    js = 0
    # include the number of genomes for each participating pvog
    mean_ani = 0
    mean_aai = 0
    
    pvogA = interaction_tuple[0]
    pvogB = interaction_tuple[1]
    
    genomesA = set((hmm_results[pvogA].keys()))
    genomesB = set(hmm_results[pvogB].keys())
    
    common_genomes = genomesA.intersection(genomesB)
    if len(common_genomes) > 0:
        all_genomes = genomesA.union(genomesB)
        
        # Jaccard score
        js = len(common_genomes) / len(all_genomes)

        # Mean ANI
        mean_ani = get_mean_from_df(common_genomes, ani_df)
        
        # Mean AAI
        mean_aai = get_mean_from_df(common_genomes, aai_df)

        # Distances
        sum_of_distances = 0
    
        for genome in common_genomes:
            hitA = hmm_results[pvogA][genome]
            hitB = hmm_results[pvogB][genome]           
            
            # Get the proper starts for distance calculation
            if hitA[-1] == '+':
                startA = hitA[1]
                endA = hitA[2]
            else: 
                startA = hitA[2]
                endA = hitA[1]
                
            if hitB[-1] == '+':
                startB = hitB[1]
                endB = hitB[2]
            else:
                startB = hitB[2]
                endB = hitB[1]            
            ## 1. Start to start
            #sum_of_distances += abs(startA - startB)
            
            ## 2. Shortest distance
            sum_of_distances += get_shortest_distance(startA, endA, startB, endB)
            
            # If they have the same orientation
            # Regardless of '+' or '-'
            if hitA[-1] == hitB[-1]:
                same_score += 1 
                
            if hitA[-1] != hitB[-1]:
                dstarts = abs(startA - startB)
                dends = abs(endA - endB)
                if dstarts >= dends:
                    inwards_score += 1                    
                else:
                    outwards_score += 1
    
        same_score = same_score / len(common_genomes)
        inwards_score = inwards_score / len(common_genomes)
        outwards_score = outwards_score / len(common_genomes)
        avg_distance = sum_of_distances / len(common_genomes)
        
    return len(genomesA), len(genomesB), len(common_genomes), js, same_score, inwards_score, outwards_score, avg_distance, mean_ani, mean_aai

#def get_ani_for_genomes(genomes, ani_df):
#    """
#    Calculate mean ani for the genomes list from the 
#    ani df
#    """
#    genomes_df = ani_df.loc[genomes, genomes]
#    return genomes_df.values.mean()

#def get_aai_for_genomes(genomes, aai_df):
#    genomes_df = aai_df.loc[genomes, genomes]
#    return genomes_df.values.mean()

def main():
    # Read in the arguments
    args = parse_args()
    print("Loading data...")

    # Store the genome sizes
    genome_sizes = get_seq_sizes(args.genomes_fasta)
    print("Loaded sequence size info for {} input sequences".format(len(genome_sizes)))
    
    print("Reading hmmsearch information...")
    # Get the hmmsearch results info for all pvogs in
    hmm_results = collect_hmmsearch_info(args.hmmer_in, genome_sizes)
    print("Done!")

    print("Reading ANI matrix...")
    # Read in the matrices
    ani_df = pd.read_csv(args.ani_matrix, index_col=0, header=0, sep="\t")
    print("Done!")

    print("Reading AAI matrix...")
    aai_df = pd.read_csv(args.aai_matrix, index_col=0, header=0, sep="\t")
    print("Done!")

    # Create a list that holds all pvog ids
    all_pvogs = get_pvogs_ids(args.profiles_file)
    
    # Create all pairs of pvogs
    all_combos = list(combinations(all_pvogs, 2))
    print("Created {} possible combinations".format(len(all_combos)))

    # Calculate scores
    print("Caclulating and writing to file...")

    counter = 0
    with open(args.outfile, 'w') as fout:
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('pvog1', 'pvog2',
                                                                   'genomes1', 'genomes2',
                                                                   'overlap_genomes', 'jaccard_score',
                                                                    'same_score', 'inwards_score', 'outwards_score',
                                                                   'avg_distance', 'mean_ani', 'mean_aai'))
        for combo in all_combos:
            int_string = '{}\t{}\t'.format(combo[0], combo[1])
            scores = calculate_scores(combo, hmm_results, ani_df, aai_df)
            scores_string = '\t'.join(map(str, scores))
            fout.write(int_string + scores_string + '\n')
            counter += 1
            if counter % 1000000 == 0:
                print("{}/{} processed".format(counter, len(all_combos)))

    print("Scores are written to {}".format(str(args.outfile)))


if __name__ == '__main__':
    main()

