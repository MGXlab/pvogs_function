#!/usr/bin/env python

import argparse
import pathlib
import gzip
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
               description="Split a multi-fasta to single files per sequence"
    )
    optionalArgs = parser._action_groups.pop()

    requiredArgs = parser.add_argument_group("required arguments")

    requiredArgs.add_argument(
        "-i",
        "--input",
        dest="input_fp",
        required=True,
        help="Input fasta. Can be gz",
    )
    requiredArgs.add_argument(
        "-o",
        "--outdir",
        dest="out_dir",
        required=True,
        help="A directory to store the files. It is NOT created",
    )
    optionalArgs.add_argument(
        "--write-reflist",
        action="store_true",
        required=False,
        dest="write_reflist",
        help="Write a file that contains all paths to the output fasta files, "
        "one per line, in the parent directory "
    )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()

def is_gz(path_string):
    """
    Return true if gzipped file
    :param path: path to file
    :return: boolean
    """
    return path_string.endswith(".gz") or path_string.endswith(".z")


def optionally_compressed_handle(path, mode):
    """
    Return a file handle that is optionally gzip compressed
    :param path: path
    :param mode: mode
    :return: handle
    """
    if mode == "r" or mode == "rb":
        mode = "rt"
    if mode == "w" or mode == "wb":
        mode = "wt"
    if is_gz(path):
        return gzip.open(path, mode=mode)
    else:
        return open(path, mode=mode)

def split_multifasta(input_fp, output_dir, write_reflist=False):
    record_no = 0
    filenames = []
    with optionally_compressed_handle(str(input_fp), 'r') as fin:
        for record in SeqIO.parse(fin, "fasta"):
            genome_acc = record.id
            single_fasta = "{}.fasta".format(genome_acc)
            single_fasta_fp = output_dir.joinpath(single_fasta)
            with open(single_fasta_fp, 'w') as fout:
                record_no += SeqIO.write(record, fout, "fasta")
            filenames.append(single_fasta_fp)

            if record_no % 10000 == 0:
                print("processed {} records".format(record_no))
    if write_reflist:
        reflist_txt = output_dir.parent.joinpath("reflist.txt")
        with open(reflist_txt, 'w') as refout:
            for f in filenames:
                refout.write('{}\n'.format(f))

    return record_no

def main():
    args = parse_args()
    a = split_multifasta(pathlib.Path(args.input_fp),
                         pathlib.Path(args.out_dir),
                         args.write_reflist)


if __name__ == '__main__':
    main()
