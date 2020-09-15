#!/usr/bin/env python

import argparse

# This is to handle ete3 plotting with Qt in a headless environment
import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
from ete3 import NCBITaxa, TextFace, TreeStyle, faces


from collections import namedtuple
import pandas as pd
from Bio import SwissProt
from pathlib import Path


parser = argparse.ArgumentParser(description="A script that parses IntAct.txt (filtered for viruses) "
                                             "to produce a visual summary of the taxa represented in it\n",
                                formatter_class=argparse.RawTextHelpFormatter)
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-o', '--output_dir',
                    dest='output_dir',
                    type=str,
                    required=True,
                    help="The output directory where all resulting files will be stored. "
                         "It is created if it's not there")
requiredArgs.add_argument('-i', '--intact_in',
                    dest='intact_in',
                    help="A txt file containing the filtered IntAct data for viruses",
                    required=True)
optionalArgs.add_argument("--tax-db",
                          dest="tax_db",
                          required=False,
                          help="Path to taxa.sqlite for ete3"
                          )

parser._action_groups.append(optionalArgs)


## Functions for parsing IntAct fields

def parse_field(field_string):
    """
    Args:
        field_string (str) : A string of the form 'xref:value(desc)'.
                            Desc is optional
    """
    field_data = {'xref': None}

    parts = field_string.split(":")
    xref = parts[0]
    field_data['xref'] = {xref: {'value': None, 'desc': None}}
    try:
        rest = parts[1]
        if "(" in rest:
            subfields = rest.split("(")
            field_data['xref'][xref]['value'] = subfields[0]
            field_data['xref'][xref]['desc'] = subfields[1][:-1]
        else:
            field_data['xref'][xref]['value'] = rest
            field_data['xref'][xref]['desc'] = '-'
    except:
        print(field_data)
    return field_data


def parse_column_string(astring):
    """
    Args:
        astring (str): A full string after splitting the PSI-MITAB on tabs
    Returns:
        column_data (dict): A dictionary of the form {'xref':
                                                        {xref:
                                                        {'value': value,
                                                        'desc': desc}
                                                        }
                                                     }
    """
    column_data = {}
    if "|" in astring:
        fields = astring.split("|")
    else:
        fields = [astring]  # LAAAAZY

    for f in fields:
        result = parse_field(f)
        k = list(result['xref'].keys())[0]
        if 'xref' not in column_data:
            column_data['xref'] = {k: result['xref'][k]}
        else:
            column_data['xref'][k] = result['xref'][k]
    return column_data


def get_single_value_from_xref_dic(xref_dic):
    """
    If the given dictionary contains only one entry, get its value,
    independent of key.

    Args:
        xref_dic (dict): A nested dictionary of the form {'xref':
                                                            {xref:
                                                            {'value': value,
                                                             'desc': desc}
                                                            }
                                                         }
    Returns:
        value (str): The value of 'value' the innermost dict
    """
    if xref_dic:
        if len(xref_dic['xref'].keys()) > 1:
            return None
        else:
            k = next(iter(xref_dic['xref']))
            value = xref_dic['xref'][k]['value']
    else:
        value = '-'
    return value


def parse_intact_data(int_file, skip_header=True):
    ProtInfo = namedtuple('ProtInfo', ['id', 'srcdb', 'taxid'])
    interactions_data = {}

    with open(int_file, 'r') as f:
        if skip_header:
            next(f)
        for line in f:
            # Make some parsable fields
            columns = [column.strip() for column in line.split('\t')]

            # Get the proteinA info
            protA_data = parse_column_string(columns[0])

            protA_id = get_single_value_from_xref_dic(protA_data)
            sourceA_db = list(protA_data['xref'].keys())[0]
            taxA_data = parse_column_string(columns[9])
            taxidA = int(taxA_data['xref']['taxid']['value'])
            protAraw = ProtInfo(id=protA_id, srcdb=sourceA_db, taxid=taxidA)

            # Get the proteinB info
            protB_data = parse_column_string(columns[1])
            protB_id = get_single_value_from_xref_dic(protB_data)
            sourceB_db = list(protB_data['xref'].keys())[0]
            taxB_data = parse_column_string(columns[10])
            taxidB = int(taxB_data['xref']['taxid']['value'])
            protBraw = ProtInfo(id=protB_id, srcdb=sourceB_db, taxid=taxidB)


            # Sorting on interactors to make unique keys for the same interaction
            # otherwise (A,B) != (B,A)
            interaction = tuple(sorted((protA_id, protB_id,)))

            if interaction[0] == protAraw.id:
                protAinfo, protBinfo = protAraw, protBraw
            else:
                protAinfo, protBinfo = protBraw, protAraw

            if protAinfo.taxid == protBinfo.taxid:
                same_taxid = 1
            else:
                same_taxid = 0

            if protAinfo.id == protBinfo.id:
                same_protein = 1
            else:
                same_protein = 0

            if str(interaction) not in interactions_data:
                interactions_data[str(interaction)] = {'no_of_evidence': 1,
                                                   'prot_A': protAinfo.id,
                                                   'prot_B': protBinfo.id,
                                                   'source_A': protAinfo.srcdb,
                                                   'source_B': protBinfo.srcdb,
                                                   'taxid_A': protAinfo.taxid,
                                                   'taxid_B': protBinfo.taxid,
                                                   'same_taxid': same_taxid,
                                                   'same_protein': same_protein}

            elif str(interaction) in interactions_data:
                interactions_data[str(interaction)]['no_of_evidence'] += 1

    return interactions_data





def write_response_text_to_file(response_fn, response_text):
    with open(response_fn, 'w') as fh:
        fh.write(response_text)


def sync_query_list_with_response(response_fn, query_list):
    db_data = {}
    with open(response_fn, 'r') as fh:
        for record in SwissProt.parse(fh):
            acc = record.accessions[0]
            # Select only EMBL and RefSeq crossrefs
            refseq_refs, embl_refs = [], []
            for db_ref in record.cross_references:
                if db_ref[0] == 'RefSeq':
                    refseq_refs.append(db_ref[1:])
                elif db_ref[0] == 'EMBL':
                    embl_refs.append(db_ref[1:])
            db_data[acc] = {'RefSeq': refseq_refs,
                            'EMBL': embl_refs}

    # This is to handle isoforms
    # E.g. P03692 and P03692-1 can both be included in the query list
    # P03705-2 can be in the query list but not P03705
    for prot in query_list:  # For each of the original queries
        if (prot not in db_data) and ('-' in prot):  # If the query is not returned
            base_name = prot.split('-')[0]  # Search for the fist part
            if base_name in db_data:  # If it is present in the db_data
                if base_name not in query_list:  # If it's not in the original query list
                    db_data[prot] = db_data[base_name]  # Fill in the information for the corresponding protein
                    db_data.pop(base_name)  # AND remove the original part
                elif base_name in query_list:  # If the first part is in the original query
                    db_data[prot] = db_data[base_name]  # Fill in the information WITHOUT removing the original part
        elif prot not in db_data:
            print("I don't know what to do with this id: {}".format(prot))
            pass

    return db_data


def uniprots_list_to_query(uniprots_list):
    query_data = {}
    for i in uniprots_list:
        primary_acc = i.split('-')[0]
        query_data[i] = primary_acc
    return query_data


def iso_is_present(k, values):
    ind = None
    for i, v in enumerate(values):
        for m in v:
            if '[{}]'.format(k) in m:
                ind = int(i)
    return ind


def select_embl_info(k, values):
    ind = None
    for i, v in enumerate(values):
        if 'ALT' in v[2]:
            pass
        else:
            ind = int(i)
    return ind


def get_refseq_info(prot_id, refs):
    if len(refs) == 1:
        refseq_genome = refs[0][1].split()[0].strip('.')
        refseq_protein = refs[0][0]
    elif ('-' in k) and (len(refs) > 1):
        isoform_index = iso_is_present(prot_id, refs)
        if isoform_index is not None:
            refseq_genome = refs[isoform_index][1].split()[0].strip('.')
            refseq_protein = refs[isoform_index][0]
        else:
            refseq_genome = refs[0][1].split()[0].strip('.')
            refseq_protein = refs[0][0]
    else:
        refseq_genome = refs[0][1].split()[0].strip('.')
        refseq_protein = refs[0][0]
    return refseq_genome, refseq_protein


def get_embl_info(prot_id, refs):
    if len(refs) == 1:
        embl_genome = refs[0][0]
        embl_protein = refs[0][1]
    elif len(refs) > 1:
        info_index = select_embl_info(prot_id, refs)
        if info_index is not None:
            embl_genome = refs[info_index][0]
            embl_protein = refs[info_index][1]
        else:
            embl_genome = refs[0][0]
            embl_protein = refs[0][1]
    else:
        embl_genome = refs[0][0]
        embl_protein = refs[0][1]
    return embl_genome, embl_protein


def map_proteins_to_genomic_accessions(db_data):
    accessions_map = {}
    for k in db_data:
        if db_data[k]['RefSeq']:
            refs = db_data[k]['RefSeq']
            (genome_id, protein_id) = get_refseq_info(k, refs)
            accessions_map[k] = {'genome_id': genome_id,
                                 'protein_id': protein_id}
        elif not db_data[k]['RefSeq'] and db_data[k]['EMBL']:
            refs = db_data[k]['EMBL']
            (genome_id, protein_id) = get_embl_info(k, refs)
            accessions_map[k] = {'genome_id': genome_id,
                                 'protein_id': protein_id}
    return accessions_map

# ETE plotting options and functions
def layout(node):
    faces.add_face_to_node(TextFace(node.sci_name), node, 0)


def plot_taxids(taxids_list, tree_png, tree_nw, tax_db=None):
    if tax_db is not None:
        ncbi = NCBITaxa(dbfile=tax_db)
    else:
        ncbi=NCBITaxa()

    tree = ncbi.get_topology(taxids_list)
    ts = TreeStyle()
    ncbi.annotate_tree(tree, taxid_attr="sci_name")
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.layout_fn = layout
    tree.render(tree_png, tree_style=ts)
    tree.write(format=1, outfile=tree_nw)


def write_sequence_dic_to_file(accession_map, fileout):
    header_string = '\t'.join(['uniprot_id', 'genome_id', 'protein_id'])
    with open(fileout, 'w') as f:
        f.write(header_string+'\n')
        for acc in accession_map:
            info_string = '\t'.join([acc, accession_map[acc]['genome_id'], accession_map[acc]['protein_id']])
            f.write(info_string + '\n')


if __name__ == '__main__':
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    else:
        print("Output directory {} already exists! Contents will be overwritten".
              format(args.output_dir))

    plot_dir = output_dir.joinpath('plots')
    plot_dir.mkdir(exist_ok=True)
    tree_png = plot_dir.joinpath('taxa_tree.png')
    tree_nw = plot_dir.joinpath('taxa_tree.nw')

    # Get the intact file
    intact_raw = args.intact_in
    intact_data = parse_intact_data(intact_raw)

    # Calculate the number of initial entries
    no_of_entries = 0
    for k in intact_data.keys():
        no_of_entries += intact_data[k]['no_of_evidence']

    print("No. of entries in input: {}".format(no_of_entries))
    print("No. of interactions: {}".format(len(intact_data.keys())))

    # Put the data in a data frame
    df = pd.DataFrame.from_dict(intact_data, orient='index')
    df.index.name = 'interaction'
    metadata_tsv = output_dir.joinpath('metadata.tsv')
    df.to_csv(metadata_tsv, sep='\t')

    # Interactors information
    unique_protAs = set(df['prot_A'])
    unique_protBs = set(df['prot_B'])
    unique_prots = len(set.union(unique_protAs, unique_protBs))
    print("No. of unique proteins: {}".format(unique_prots))

    # Taxids information
    unique_taxidsA = set(df['taxid_A'])
    unique_taxidsB = set(df['taxid_B'])
    unique_taxids = set.union(unique_taxidsA, unique_taxidsB)
    print("No. of unique taxids: {}".format(len(unique_taxids)))

    print("Plotting NCBI taxonomy information. Results are stored in: {}".format(plot_dir))
    plot_taxids(unique_taxids, tree_png.as_posix(), tree_nw.as_posix(), args.tax_db)

    ######################################
    ### UP TO HERE EVERYTHING IS RAW
    ### NOW WE START SUBSETTING
    #######################################
    # Getting the uniprot related entries to dispatch a query
    # protAs = df.loc[(df['source_A'] == 'uniprotkb')]['prot_A'].unique().tolist()
    # protBs = df.loc[(df['source_B'] == 'uniprotkb')]['prot_B'].unique().tolist()
    # uniprots = set.union(set(protAs), set(protBs))
    # uniprots_list = list(uniprots)
    #
    # print("{} identifiers will be searched against uniprot for genome and protein sequence retrieval".
    #       format(len(uniprots_list)))
    #
    # uniprot_response = output_dir.joinpath('uniprot_info.txt')
    # print("Dispatching queries...")
    # final_result = []
    # # Split the query into batches of 100
    # for i in range(0, len(uniprots_list), 100):
    #     slise = uniprots_list[i:i + 100]
    #     # Run the query for each batch
    #     result = map_retrieve(slise)
    #     # And store the result in a list
    #     final_result.append(result)
    #
    # # Store the result in a file
    # response_text = ''.join(final_result)
    # write_response_text_to_file(uniprot_response, response_text)
    #
    # db_data = sync_query_list_with_response(uniprot_response, uniprots_list)
    # # For each original entry in the protein list,
    # # store its primary accession
    # # TO DO
    # # the sync_query_list_with_response function needs revisiting
    # query_data = uniprots_list_to_query(uniprots_list)
    #
    # accession_map = map_proteins_to_genomic_accessions(db_data)
    # for protein in query_data:  # {alt_acc : primary_acc}
    #     if protein not in accession_map:  # alt_acc not in results
    #         try:
    #             primary_acc = query_data[protein]  # get the primary accession
    #             accession_map[protein] = accession_map[primary_acc]
    #         except KeyError:
    #             accession_map[protein] = {'genome_id': '-', 'protein_id': '-'}
    #
    # sequence_info = output_dir.joinpath('sequence_info.txt')
    # write_sequence_dic_to_file(accession_map, sequence_info)
