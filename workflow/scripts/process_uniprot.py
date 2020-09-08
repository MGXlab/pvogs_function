#!/usr/bin/env python

import argparse
from Bio import SwissProt, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pathlib
from ast import literal_eval as make_tuple
import pandas as pd

parser = argparse.ArgumentParser(description="This is to get several mapping files for the interactions.",
                                 formatter_class=argparse.RawTextHelpFormatter)
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--interactions-tsv',
                          dest='interactions_tsv',
                          required=True,
                          help="Filtered interactions_tsv from get_uniprot_ids.py")
requiredArgs.add_argument('-l', '--input_list',
                          dest='uniprots_list',
                          type=str,
                          required=True,
                          help="The output of get_uniprot_ids.py")
requiredArgs.add_argument('-s', '--swissprot-file',
                          dest='swissprot_file',
                          required=True,
                          help="A SwissProt file, containing results based on the primary ids")
requiredArgs.add_argument('-p', '--prefix-out',
                          dest='prefix_out',
                          required=True,
                          help="The prefix of the file names to write to. This will produce files:\n"
                               "<prefix>/proteins.faa,\n"
                               "<prefix>/uniprot2ncbi.mapping.txt\n"
                               "<prefix>/ncbi2uniprot.mapping.txt\n"
                               "<prefix>/skipped.uniprot.txt\n"
                               "<prefix>/ncbi_interactions.txt\n"
                               "<prefix>/ncbi_genomes.txt")
parser._action_groups.append(optionalArgs)


def has_refseq(db_list):
    """
    Return the index of the list where the 'RefSeq' string is located.
    Otherwise return None
    :param db_list: A list of db names taken as the first element of the tuples in a
    Swissprot.record.cross_references list
    :return: int: index or None
    """
    if 'RefSeq' in db_list:
        return db_list.index('RefSeq')
    else:
        return None


def give_me_proper_embl(cross_refs):
    """
    Filter for references where the first element == 'EMBL',
    then search for the first occurence where the genome accession is not '-'.
    This is to get both a valid protein accession and genome accession.

    :param cross_refs: The full list of SwissProt.record.cross_references
    :return:
    """
    # Get embl records first
    embl_refs = []
    for ref in cross_refs:
        if ref[0] == 'EMBL':
            embl_refs.append(ref)

    genome_acc, prot_acc = embl_refs[0][1], embl_refs[0][2]

    for ref in embl_refs:
        if ref[2] != '-':
            genome_acc, prot_acc = ref[1], ref[2]
            break
    return genome_acc, prot_acc


def has_embl(db_list):
    """
    Check if 'EMBL' is in a list of database names
    :param db_list:
    :return: True or None
    """
    if 'EMBL' in db_list:
        return True
    else:
        return None


def get_ncbi_accessions(cross_refs):
    """
    Given a SwissProt.record.cross_references list, extract a tuple
    of genome,protein accessions
    :param cross_refs: Bio.SwissProt.record.cross_references list
    :return: (genome_accession, protein_accession) associated with the record
    """
    dbs = [ref[0] for ref in cross_refs]

    if has_refseq(dbs) is not None:
        loc = has_refseq(dbs)
        genome_acc = cross_refs[loc][2]
        if '[' in genome_acc:
            genome_acc = genome_acc.split()[0].rstrip('.')
        prot_acc = cross_refs[loc][1]
    elif has_embl(dbs) is not None and has_refseq(dbs) is None:
        genome_acc, prot_acc = give_me_proper_embl(cross_refs)
    #         genome_acc, prot_acc = cross_refs[loc][1], cross_refs[loc][2]
    else:
        genome_acc, prot_acc = '-', '-'

    return genome_acc, prot_acc


def construct_description(protein_seq, uniprotkb_id, ncbi_genome, ncbi_protein, start=None, end=None, taxid=None):
    if start is not None and end is not None:
        return f'uniprotkb={uniprotkb_id};' \
            f'ncbi_protein={ncbi_protein};' \
            f'start={start};' \
            f'end={end};' \
            f'ncbi_genome={ncbi_genome};' \
            f'taxid={taxid}'
    else:
        return f'uniprotkb={uniprotkb_id};' \
            f'ncbi_protein={ncbi_protein};' \
            f'start=1;' \
            f'end={len(protein_seq)};' \
            f'ncbi_genome={ncbi_genome};' \
            f'taxid={taxid};'


def extract_prochain_sequence_from_record(record,
                                          primary_id,
                                          prochain_id,
                                          genome_acc='-',
                                          protein_acc='-'):
    r_taxid = record.taxonomy_id
    for feature in record.features:
        if feature[0] == 'CHAIN' and feature[-1] == prochain_id:
            start, end = feature[1], feature[2]
            rec_seq = Seq(record.sequence[start - 1:end], IUPAC.protein)
            original = '-'.join([primary_id, prochain_id])
            if protein_acc != '-':
                new_id = f'{protein_acc}_at_{start}-{end}'
                prot_rec = SeqRecord(rec_seq, id=new_id,
                                     description=construct_description(str(rec_seq),
                                                                       original,
                                                                       genome_acc,
                                                                       protein_acc,
                                                                       start,
                                                                       end,
                                                                       taxid=r_taxid))
            else:
                prot_rec = SeqRecord(rec_seq, id=original,
                                     description=construct_description(str(rec_seq),
                                                                       original,
                                                                       genome_acc,
                                                                       protein_acc,
                                                                       start,
                                                                       end,
                                                                       taxid=r_taxid))
    return prot_rec

def uniprot_id_type(uniprot_id):
    """
    Check the type of a given uniprot_id
    P04591 -> 'primary'
    P04591-PRO_0000038593 -> 'prochain'
    P04591-3 -> 'isoform
    :param uniprot_id: String of uniprot_id
    :return:
    """
    if '-' in uniprot_id:
        second_part = uniprot_id.split('-')[1]
        if second_part.startswith('PRO'):
            return 'prochain'
        else:
            return 'isoform'
    else:
        return 'primary'

def give_me_the_record(primary_id, swissprot_file):
    """
    Return a single record given with the primary id
    :param primary_id: A primary id
    :param swissprot_file: A swissprot file
    :return: A record with accession == primary id
    """
    with open(swissprot_file, 'r') as fh:
        for record in SwissProt.parse(fh):
            if primary_id in record.accessions:
                return record

def extract_protein_from_record(record):
    """
    Grab the protein sequence as a string from a SwissProt record
    :param record: A Bio.SwissProt.SeqRecord instance
    :return:
    """
    return str(record.sequence)


if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.uniprots_list, 'r') as fin:
        uniprots_list = [line.strip() for line in fin]

    # Create the names of the output files
    output_dir = pathlib.Path(args.prefix_out)
    output_dir.mkdir(exist_ok=True)
    # This is the fasta file that contains the proteins
    output_fasta = pathlib.Path.joinpath(output_dir, 'proteins.faa')


    # Initialize several book-keeping lists and dicts
    # Some records I skip
    skipped = []
    # A dict {uniprotkb : ncbi-accession , ...}
    uniprot2ncbi_protein_mapping = {}
    # Keep track of how many sequences are being written
    written = 0
    # A dict { ncbi_protein_accession : [ uniprot_accession1, ...] , ... }
    duplicates = {}
    genomes = []

    with open(output_fasta, 'w') as fout:
        # Loop over the initial list of interactors
        for protein_id in uniprots_list:
            # Is it primary/prochain/isoform?
            entry_type= uniprot_id_type(protein_id)
            if entry_type == 'primary':
                # Grab the associated record
                r = give_me_the_record(protein_id, args.swissprot_file)
                # Get the genome an protein accessions
                ncbi_genome, ncbi_protein = get_ncbi_accessions(r.cross_references)

                # I want to map the uniprot ids to ncbi ones.
                # If there is nothing to map skip the record altogether
                if ncbi_protein == '-':
                    skipped.append(protein_id)
                    pass
                else:
                    # For primary accessions the protein is the same
                    rec_id = ncbi_protein
                    rec_seq = extract_protein_from_record(r)
                    rec_taxid = r.taxonomy_id
                    description = construct_description(rec_seq, protein_id, ncbi_genome, ncbi_protein, taxid=rec_taxid)
                    # Construct a SeqRecord object to write to file
                    prot_rec = SeqRecord(Seq(rec_seq, IUPAC.protein),
                                         id=rec_id,
                                         description=description)
                    # Book-keeping
                    # match uniprot id to its ncbi counterpart
                    uniprot2ncbi_protein_mapping[protein_id] = prot_rec.id
                    # get the genome
                    genomes.append(ncbi_genome)

                    # Avoid duplicates
                    if prot_rec.id not in duplicates:
                        # Put this unique identifier in the seen_proteins dict
                        duplicates[prot_rec.id] = [protein_id]
                        # And write to file
                        written += SeqIO.write(prot_rec, fout, "fasta")
                    else:
                        print("Possible duplicate:{}".format(prot_rec.id))
                        duplicates[prot_rec.id].append(protein_id)
                        pass
            elif entry_type == 'prochain':
                # When a PRO is present, things need to be different
                primary_id = protein_id.split('-')[0]
                prochain = protein_id.split('-')[1]
                r = give_me_the_record(primary_id, args.swissprot_file)
                ncbi_genome, ncbi_protein = get_ncbi_accessions(r.cross_references)
                if ncbi_protein == '-':
                    skipped.append(protein_id)
                    pass
                else:
                    # The actual sequence is in the features table
                    prot_rec = extract_prochain_sequence_from_record(r, primary_id, prochain,
                                                                     genome_acc=ncbi_genome,
                                                                     protein_acc=ncbi_protein)
                    uniprot2ncbi_protein_mapping[protein_id] = prot_rec.id
                    genomes.append(ncbi_genome)

                    # Avoid duplicates
                    if prot_rec.id not in duplicates:
                        written += SeqIO.write(prot_rec, fout, "fasta")
                        duplicates[prot_rec.id] = [protein_id]
                    else:
                        print("Possible duplicate {}".format(prot_rec.id))
                        duplicates[prot_rec.id].append(protein_id)
                        pass
            elif entry_type == 'isoform':
                skipped.append(protein_id)
                print("{} is an isoform and I do not handle these for now".format(protein_id))
                pass
    print(30*'=')
    print("{} proteins were extracted".format(written))

    # This is a file that contains for each uniprot id its mapping to an ncbi protein
    uniprot2ncbi_mapping_txt = pathlib.Path.joinpath(output_dir, 'uniprot2ncbi.mapping.txt')
    uniprot_counter = 0
    with open(uniprot2ncbi_mapping_txt, 'w') as fout:
        for uniprot_id, ncbi_id in uniprot2ncbi_protein_mapping.items():
            uniprot_counter += 1
            fout.write('{}\t{}\n'.format(uniprot_id, ncbi_id))
    print('{} uniprot ids are mapped to NCBI accessions'.format(uniprot_counter))

    # A list of genomes to download from NCBI
    genome_accessions_txt = pathlib.Path.joinpath(output_dir, 'genome_accessions.txt')
    with open(genome_accessions_txt, 'w') as fout:
        for genome in set(genomes):
            fout.write('{}\n'.format(genome))
    print('{} genomes are written to file'.format(len(set(genomes))))

    # Each ncbi protein, might map to more uniprot ids
    # This files contains this information - an aggregation at ncbi level
    ncbi2uniprot_mapping_txt = pathlib.Path.joinpath(output_dir, 'ncbi2uniprot.mapping.txt')
    with open(ncbi2uniprot_mapping_txt, 'w') as fout:
        for ncbi_p in duplicates:
            fout.write('{}\t{}\n'.format(ncbi_p, ','.join(duplicates.get(ncbi_p))))

    # This files contains only the duplicates
    duplicates_tsv = pathlib.Path.joinpath(output_dir, 'duplicates.tsv')
    with open(duplicates_tsv, 'w') as fout:
        for ncbi_p in duplicates:
            if len(duplicates.get(ncbi_p)) > 1:
                fout.write('{}\t{}\n'.format(ncbi_p, ','.join(duplicates.get(ncbi_p))))

    # This files contains the list of records that were originally input
    # but were skipped. For now this happens
    # (1) when no valid ncbi_protein accession is found
    # (2) for isoforms
    skipped_txt = pathlib.Path.joinpath(output_dir, 'skipped.txt')
    with open(skipped_txt, 'w') as fout:
        for uniprot in skipped:
            fout.write('{}\n'.format(uniprot))

    # Get the proteins that are actually covered through this filtering process
    proteins_covered = []
    for ncbi_p in duplicates:
        for uniprot_p in duplicates.get(ncbi_p):
            proteins_covered.append(uniprot_p)
    print('{} uniprot ids are covered'.format(len(proteins_covered)))

    # Get the interactions from the input file
    interactions = pd.read_csv(args.interactions_tsv, sep='\t')
    print('Input file contains {} interactions'.format(interactions.shape[0]))
    # Filter the interactions dataframe to these interactions only
    interactions_filtered = interactions.loc[interactions.prot_A.isin(proteins_covered)
                                             & interactions.prot_B.isin(proteins_covered), ]
    print('{} interactions are covered from the {} proteins above'.format(interactions_filtered.shape[0],
                                                                          len(proteins_covered)))

    # Write the final set of interactions to file
    interactions_filtered_tsv = pathlib.Path.joinpath(output_dir, 'interactions_filtered.tsv')
    interactions_filtered.to_csv(interactions_filtered_tsv, sep='\t', index=False)

    # Also provide a tsv file that contains the uniprot-uniprot interaction
    # as ncbi-ncbi interaction (for the negatives set construction)
    ncbi_interactions_tsv = pathlib.Path.joinpath(output_dir, 'ncbi_interactions.tsv')
    # This is a list of strings ["('intA', 'intB')" , ...]
    interaction_pairs = interactions_filtered['interaction'].to_list()
    with open(ncbi_interactions_tsv, 'w') as fout:
        for pair in interaction_pairs:
            # Convert the literal tuple string to actual tuple
            # TO DO
            # Just make it a tuple to begin with in the interactions file...
            p = make_tuple(pair)
            fout.write('{}\t{}\t{}\t{}\n'.format(p[0], p[1],
                                                 uniprot2ncbi_protein_mapping.get(p[0]),
                                                 uniprot2ncbi_protein_mapping.get(p[1])))
    print('Done!')
