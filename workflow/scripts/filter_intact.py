#!/usr/bin/env python
from ete3 import NCBITaxa, TextFace, TreeStyle, faces
import argparse


parser = argparse.ArgumentParser(description='Parse the IntAct DB text file to get viral related entries')

# Store the default optional groups to the optionalArgs
optionalArgs = parser._action_groups.pop()

# Add requiredArgs as an argument group for displaying help better
requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', dest="intact_in",
                    help="The raw txt intact file",
                    required=True)
requiredArgs.add_argument('-o', dest="intact_out",
                    required=True,
                    help="File path to write results to")

# Add the rest of the optional arguments
optionalArgs.add_argument('-v', '--vhost',
                    dest="vhost",
                    required=False,
                    help="A txt file containing a list of bacteria and archaea associated virus taxids")
optionalArgs.add_argument('--phages-only',
                          dest='phages_only',
                          action='store_true',
                         help='Filter taxids to only contain phages provided in the vhost file')
optionalArgs.add_argument("--tax-db",
                          dest="tax_db",
                          required=False,
                          help="Location of the NCBI taxonomy sqlite db, created"
                          )
# Add the optionalArgs to the action_groups
parser._action_groups.append(optionalArgs)

def file_to_list(filename):
    """
    Read in a one-column txt file to a list
    :param filename:
    :return: A list where each line is an element
    """
    with open(filename, 'r') as fin:
        alist = [line.strip() for line in fin]
    return alist


def parse_field(field_string):
    """
    Get a tuple of the first entries found in an psimitab column
    (xref, value, desc)
    :param field_string:
    :return: A tuple (xref, value, desc) - desc is optional
    """
    parts = field_string.split(":")
    xref = parts[0]
    rest = parts[1]
    if "(" in rest:
        subfields = rest.split("(")
        value = subfields[0]
        desc = subfields[1][:-1]
    else:
        value = rest
        desc = "NA"
    return xref, value, desc


def parse_column_string(astring):
    """
    Parses the information stored in a column field of psi-mitab
    :param astring: The string of psimitab
                    e.g. `taxid:9606(human)|taxid:9606(Homo sapiens)`
    :return: A dictionary of the values
            {'xref': ['taxid'],
            'value': ['9606'],
            'desc': ['human', 'Homo sapiens']}
    """

    # Initialize an empty list for each xref, value, desc key
    column_data = {k: [] for k in ['xref', 'value', 'desc']}
    if "|" in astring:
        fields = astring.split("|")
    else:
        fields=[astring]

    for f in fields:
        result = parse_field(f)
        xref, value, desc = result[0], result[1], result[2]
        if xref not in column_data['xref']:
            column_data['xref'].append(xref)
        if value not in column_data['value']:
            column_data['value'].append(value)
        if desc not in column_data['desc']:
            column_data['desc'].append(desc)
    return column_data


## TREE STUFF ##
def layout(node):
    """
    Helper function to annotate nodes with scientific name
    :param node: A node instance of a ete3.Tree
    :return: Annotated node with sci_name
    """
    faces.add_face_to_node(TextFace(node.sci_name), node, 0)


if __name__ == '__main__':
    # Parse the args
    args = parser.parse_args()
    ncbi = NCBITaxa(dbfile=args.tax_db)
    ncbi_viruses = ncbi.get_descendant_taxa('10239', intermediate_nodes=True)

    intact_in = args.intact_in
    intact_out = args.intact_out

    # If a file with host taxids is provided
    if args.vhost:
        phage_taxids_file = args.vhost
        # Get the list of taxids
        phage_taxids = list(map(int, file_to_list(phage_taxids_file)))
        print("Phages list : {} entries".format(len(phage_taxids)))
    else:
        phage_taxids = set()
    # If we want to only focus on phages
    if args.phages_only:
        ncbi_phages = ncbi.get_descendant_taxa('28883', intermediate_nodes=True) # 28883 -> Caudovirales
        compare_list = set.union(set(phage_taxids), set(ncbi_phages))
    else:
        ncbi_virus_taxids = ncbi.get_descendant_taxa('10239', intermediate_nodes=True)
        print("NCBI Viruses descendants: {}".format(len(ncbi_virus_taxids)))
        compare_list = set.union(set(ncbi_virus_taxids), set(phage_taxids))

    print("Contents of input file will be scanned against {} viral associated taxids.".
          format(len(compare_list)))
    # Initialize an interaction counter and a virus-virus counter
    int_counter, vv_counter = 0,0
    # Store the parsed taxids
    taxids = []
    with open(intact_in, 'r') as fin, open(intact_out, 'w') as fout:
        header_line = fin.readline()
        fout.write(header_line)
        for line in fin:
            int_counter += 1
            # Print some progress
            if int_counter % 10000 == 0:
                print("Parsed {} entries".format(int_counter))
            columns = [field.strip() for field in line.split("\t")]
            # On PSIMI-TAB v25 columns 8,9 have the taxonomy info
            if columns[9] != '-' and columns[10] != '-':
                taxidA_data = parse_column_string(columns[9])
                taxidB_data = parse_column_string(columns[10])
                taxidA, taxidB = taxidA_data['value'], taxidB_data['value']
                # If the taxids are both viral
                if (int(taxidA[0]) in compare_list) and (int(taxidB[0]) in compare_list):
                    vv_counter += 1
                    fout.write(line)
                    taxids.append(taxidA[0])
                    taxids.append(taxidB[0])

    print("{} entries were parsed.".format(int_counter))
    print("{} entries were viral".format(vv_counter))
    taxids_set = set(taxids)
    print("Filtered taxids list contains {} taxa".format(len(taxids_set)))
