#!/usr/bin/env python

import argparse
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from math import ceil

parser = argparse.ArgumentParser(description="Given a txt file with uniprot ids, post a query to uniprot\n"
                                             "and retrieve a SwissProt file with the result",
                                formatter_class=argparse.RawTextHelpFormatter)
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--input_txt',
                    dest='input_txt',
                    type=str,
                    required=True,
                    help="The output of get_uniprot_ids.txt")
requiredArgs.add_argument('-o', '--output_sp',
                    dest='output_sp',
                    help="A txt file containing the raw uniprot ids from intact.\n"
                         "This will contain the response from uniprot",
                    required=True)
optionalArgs.add_argument("--filter-list",
                          action="store_true",
                          help="Specify if the input list needs to be filtered for PROs and isoforms",
                          dest="filter",
                          required=False)
optionalArgs.add_argument("--output-fmt",
                          default='txt',
                          type=str,
                          dest="output_fmt",
                          help="SwissProt records are returned by default.\n"
                               "You can use 'fasta' for retrieving only the fasta sequences\n"
                               "or 'xml' to get an xml file.")

parser._action_groups.append(optionalArgs)


def get_primary_ids_from_list(uniprot_ids):
    """
    Filter the input list to return only querrable primary ids.
    As primary I define the part of the identifier's string that
    serves as the accession.
    E.g. P04591 is a primary id
        P04591-PRO_0000038593 is a PTM chain with primary id P04591
        P04591-3 is an isoform with primary id P04591

    This will only return P04591. The rest may be retrieved from parsing
    the results of the query.

    :param uniprot_ids: A list of uniprot specific identifiers
    :return: primary_ids: A list of primary ids
    """
    primary_ids = []
    for i in uniprot_ids:
        if '-' in i:
            fields = i.split('-')
            primary_id = fields[0]
            # Append the primary id to the query list
            if primary_id not in primary_ids:
                primary_ids.append(primary_id)

            # Make a mapping of primary ids to their pro-chains and isoforms
        elif i not in primary_ids:
            primary_ids.append(i)

    # Double check I didn't put any duplicates in the final list
    if len(primary_ids) != len(set(primary_ids)):
        primary_ids = list(set(primary_ids))

    return primary_ids

# SHAMELESSLY STOLEN FROM
# https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf
# and modified to default to text

def map_retrieve(ids2map, source_fmt='ACC+ID', target_fmt='ACC', output_fmt='txt'):
    """
    Map database identifiers from/to UniProt accessions.
    The mapping is achieved using the RESTful mapping service provided by UniProt.
    While a great many identifiers can be mapped the documentation
    has to be consulted to check which options there are and what the database codes are.
    Mapping UniProt to UniProt effectlvely allows batch retrieval
    of entries.
    Args:
        ids2map (list or string): identifiers to be mapped
        source_fmt (str, optional): format of identifiers to be mapped. Defaults to ACC+ID, which are UniProt accessions or IDs.
        target_fmt (str, optional): desired identifier format. Defaults to ACC, which is UniProt accessions.
        output_fmt (str, optional): return format of data. Defaults to list.
    Returns:
        mapped information (str)
    """

    BASE = 'http://www.uniprot.org'
    KB_ENDPOINT = '/uniprot/'
    TOOL_ENDPOINT = '/uploadlists/'

    if hasattr(ids2map, 'pop'):
        ids2map = ' '.join(ids2map)

    payload = {'from': source_fmt,
               'to': target_fmt,
               'format': output_fmt,
               'query': ids2map,
               'include': 'yes',  # This is to include isoforms
               }

    with requests.Session() as s:
        retries = Retry(total=3,
                        backoff_factor=0.5,
                        status_forcelist=[500, 502, 503, 504])

        s.mount('http://', HTTPAdapter(max_retries=retries))
        #         s.mount('https://', HTTPAdapter(max_retries=retries))

        response = s.get(BASE + TOOL_ENDPOINT, params=payload)

        if response.ok:
            return response.text
        else:
            print(response.url)
            response.raise_for_status()


if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.input_txt, 'r') as fin:
        uniprots_list = [line.strip() for line in fin]
    print("Input file contains {} records".format(len(uniprots_list)))

    if args.filter:
        print("Filtering out PROs and isoforms")
        uniprots_list = get_primary_ids_from_list(uniprots_list)
        print("Final list contains {} primary ids".format(len(uniprots_list)))

    print("Posting queries of 100 to uniprot")
    final_result = []
    total_batches = ceil(len(uniprots_list) / 100)
    batch_no = 0
    for i in range(0, len(uniprots_list), 100):
        batch_no += 1
        slise = uniprots_list[i:i + 100]
        # Run the query for each batch
        print("Running batch {}/{}".format(batch_no, total_batches))
        result = map_retrieve(slise, output_fmt=args.output_fmt)
        # And store the result in a list
        final_result.append(result)

    with open(args.output_sp, 'w') as fout:
        fout.write(''.join(final_result))
