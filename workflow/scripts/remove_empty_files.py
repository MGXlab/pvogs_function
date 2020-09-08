#!/usr/bin/env python

import argparse
import pathlib
import os

parser = argparse.ArgumentParser(description="Remove all files associated with a genome that are empty",
                                 formatter_class=argparse.RawTextHelpFormatter)
optionalArgs = parser._action_groups.pop()

requiredArgs = parser.add_argument_group("required arguments")
requiredArgs.add_argument('-i', '--inut-dir',
                          dest='input_dir',
                          required=True,
                          help="Input genes dir that contains the results of `comparem call_genes`")
requiredArgs.add_argument('-o', '--output-txt',
                          dest='output_txt',
                          type=str,
                          required=True,
                          help="Output txt to write result. Will contain the genomes removed, if any,\n"
                               "else it will have a 'All genomes passed!' so that is not empty")

def check_file_is_empty(fin):
    """

    :param fin: A pathlib.Path object
    :return: the file if true, nothing if false
    """
    if fin.stat().st_size == 0:
        return fin
    else:
        pass

if __name__ == '__main__':
    args = parser.parse_args()

    input_path = pathlib.Path(args.input_dir)
    remove_these = []
    for f in list(input_path.glob('*.faa')):
        if check_file_is_empty(f):
            remove_these.append(f)

    with open(args.output_txt, 'w') as fout:
        if len(remove_these) != 0:
            for f in remove_these:
                print('Removing file {}'.format(f.as_posix()))
                f.unlink()
                fout.write('{}\n'.format(f.as_posix()))
        else:
            fout.write("All genomes passed!\n")
