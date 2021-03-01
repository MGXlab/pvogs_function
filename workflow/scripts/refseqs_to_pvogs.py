#!/usr/bin/env python

from Bio import SearchIO, SeqIO
from pathlib import Path
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='From uniprot ids to pvogs')

    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")
    requiredArgs.add_argument('-i', '--interactions-file',
                              dest='int_file',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A 2-column tsv file containing "
                              " interaction data between refseq identifiers" 
                              )
    requiredArgs.add_argument('-f', '--input-fasta',
                              dest='fasta_in',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="A fasta file with protein sequences "
                              "retrieved from uniprot"
                              )
    requiredArgs.add_argument('-hmm', '--input-hmm',
                              dest='hmmer_tblout',
                              type=lambda p: Path(p).resolve(strict=True),
                              required=True,
                              help="The tblout file from hmmsearching the "
                              "proteins in <input_fasta> against all pvogs"
                              )
    requiredArgs.add_argument('-o', '--output-file',
                              dest='outfile',
                              required=True,
                              type=lambda p: Path(p).resolve(),
                              help="File path to write the results in"
                              )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()


# This is higly specific for my case...
remapper = {'P17312-2': 'P17312',
           'P17312-3': 'P17312',
           'P17312-4': 'P17312',
           'P03705': 'P03705-2',
           'P03711-2': 'P03711',
           'P03692-2': 'P03692'}


def get_interactions_info(interactions_fp):
    interactions = []
    with open(interactions_fp, 'r') as fin:
        for line in fin:
            fields = line.split('\t')
            # Specific for get_uniprot_ids.py
            intA, intB = fields[0].strip(), fields[1].strip()
            interaction = tuple(sorted([intA, intB]))
            interactions.append(interaction)
    return interactions


def get_unique_interactor_ids(interactions_list):
    interactor_ids = []
    for i in interactions_list:
        interactor_ids.extend([i[0], i[1]])
    return set(interactor_ids)


def get_ids_from_fasta_file(fasta_in):
    record_ids = []
    with open(fasta_in, 'r') as fin:
        for record in SeqIO.parse(fin, "fasta"):
            record_ids.append(record.id)
    return record_ids


def parse_hmm_table(hmm_table):
    """
    Create a dictionary that stores all results per protein
    """
    hmm_results = {}
    with open(hmm_table, 'r') as fin:
        for record in SearchIO.parse(fin, "hmmer3-tab"):
            for hit in record.hits:
                prot_acc = hit.id
                target = remapper.get(prot_acc, prot_acc)
                if target not in hmm_results:
                    hmm_results[target] = [
                      (record.id, hit.evalue, hit.bitscore)
                    ]
                else:
                    hmm_results[target].append(
                      (record.id, hit.evalue, hit.bitscore)
                    )
    return hmm_results


def translate_proteins_to_pvogs(hmm_results_dic):
    '''
    Translate refseq accessions to pvogs, based on the best scoring hit
    
    Arguments:
      hmm_results_dic: dict: A dictionary that holds all hmm results per
        protein. Of the form
          { protein_id : [
              (pvog1, evalue, bitscore),
              (pvog2, evalue, bitscore), 
              ...
            ],
            ...
          }

    Return:
      protein_to_pvog: dict: A dictionary of the form
        { protein_id : {'pvog' : pvog_id,
                        'evalue': evalue,
                        'bitscore' : bitscore
                       },
            ...
        }

    '''
    protein_to_pvog = {}
    for p in hmm_results_dic:
        results = hmm_results_dic[p]
        if len(results) == 1:
            best_result = results[0]
        else:
            # Select the best scoring pvog
            scores = [r[2] for r in results]
            best_index = scores.index(max(scores))
            best_result = results[best_index]

        pvog_id, evalue, score = best_result[0], best_result[1], best_result[2]

        protein_to_pvog[p] = {
            'pvog': pvog_id,
            'bitscore' : score,
            'evalue' : evalue,
            }
    
    return protein_to_pvog


def main():
    args = parse_args()

    interactions = get_interactions_info(args.int_file)
    print("Parsed {} interactions".format(len(interactions)))

    unique_ids = get_unique_interactor_ids(interactions)
    print("Unique ids: {}".format(len(unique_ids)))

    fasta_ids = get_ids_from_fasta_file(args.fasta_in)
    print("Entries in fasta: {}".format(len(fasta_ids)))

    hmm_results = parse_hmm_table(args.hmmer_tblout)
    print("Proteins with hmmer hits: {}".format(len(list(hmm_results.keys()))))

    protein_pvog_map = translate_proteins_to_pvogs(hmm_results)
    for p in protein_pvog_map:
        print(p, protein_pvog_map[p]['pvog'], 
                protein_pvog_map[p]['bitscore'],
                protein_pvog_map[p]['evalue']
            )
    no_results = []
    counter = 0
    with open(args.outfile, 'w') as fout:
        for inter in interactions:
            intA, intB = inter[0], inter[1]
            if intA in protein_pvog_map:
                if intB in protein_pvog_map:
                    fout.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                intA, 
                                protein_pvog_map[intA]['pvog'],
                                protein_pvog_map[intA]['evalue'],
                                protein_pvog_map[intA]['bitscore'],
                                intB, 
                                protein_pvog_map[intB]['pvog'],
                                protein_pvog_map[intB]['evalue'],
                                protein_pvog_map[intB]['bitscore'],
                                )
                            )
                    counter += 1
                else:
                    no_results.append(intB)
            else:
                no_results.append(intA)

    print("Interactions translated to pvogs: {}".format(counter))
    print("Proteins with no hits: {}".format(len(set(no_results))))

    if len(no_results) > 0:
        print(20 * '=')
        for i in no_results:
            print(i)
        print(20 * '=')

    print("Results written in {}".format(args.outfile.__str__()))


if __name__ == '__main__':
    main()
