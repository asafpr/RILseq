#!/usr/bin/env python

"""
Read the list of chimeric interactions and generate a file that can be read
by circos.
"""

import sys
import argparse
from collections import defaultdict
from math import log

import RILseq

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Generate circos data file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'reads_in',
        help='An output file of map_chimeric_fragments.py with the chimeric'
        ' fragments.')
    parser.add_argument(
        '-r', '--region', type=int, default=200,
        help='Split the genome to windows of this size.')
    parser.add_argument(
        '-c', '--chrn', default='chr',
        help='Name of chromosome to plot.')
    parser.add_argument(
        '-p', '--print_chr', default='ecmain',
        help='Name of chromosome in circos.')
    parser.add_argument(
        '-m', '--min_interactions', type=int, default=100,
        help='Minimum number of interactions between two regions to plot.')
    settings = parser.parse_args(argv)

    return settings

def main(argv=None):
    settings = process_command_line(argv)
    region_interactions, _, _, _=\
        RILseq.read_reads_table(open(settings.reads_in), settings.region)
    both_strs = defaultdict(lambda: defaultdict(int))
    for reg1 in region_interactions:
        if reg1[2] != settings.chrn:
            continue
        for reg2 in region_interactions[reg1]:
            if reg2[2] != settings.chrn:
                continue
            both_strs[reg1[0]][reg2[0]] += len(region_interactions[reg1][reg2])
    for r1 in both_strs:
        for r2 in both_strs[r1]:
            if both_strs[r1][r2] > settings.min_interactions:
                sys.stdout.write('%s %d %d %s %d %d thickness=%dp\n'%(
                        settings.print_chr, r1+1, r1+settings.region,
                        settings.print_chr, r2+1, r2+settings.region,
                        log(both_strs[r1][r2])/log(10)))

    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
