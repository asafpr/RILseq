#!/usr/bin/env python

"""
Read the list of chimeric interactions and generate a file that can be read
by circos.
"""

import sys
import argparse
from collections import defaultdict
from math import log
import gzip

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
        '-s', '--summary',
        help='Plot only significant interactions that appear in the summary file.')
    parser.add_argument(
        '-r', '--region', type=int, default=200,
        help='Split the genome to windows of this size.')
    parser.add_argument(
        '--sRNAs', default=False, action='store_true',
        help='Color the lines going to or coming from sRNAs in orange. Must be'
        ' used with ec_dir.')
    parser.add_argument(
        '--known',
        help='Use this file to color the known interactions in red.'
        ' Must give --refseq_dir as well.')
    parser.add_argument(
        '--refseq_dir', default='/home/users/assafp/EC/',
        help='RefSeq dir of organism to get the gene description from.'
        ' Should be given if --known is given.')
    parser.add_argument(
        '--ec_dir', #default='/home/users/assafp/Database/EcoCyc/19.0/data',
        help='EcoCyc data dir, used to map the regions to genes. If not'
        ' given only the regions will be reported')
    parser.add_argument(
        '-c', '--chrn', default='chr',
        help='Name of chromosome to plot.')
    parser.add_argument(
        '-p', '--print_chr', default='ecmain',
        help='Name of chromosome in circos.')
    parser.add_argument(
        '-m', '--min_interactions', type=int, default=100,
        help='Minimum number of interactions between two regions to plot.')
    parser.add_argument(
        '--EC_chrlist', default='COLI-K12,chr',
        help='A comma separated dictionary of chromosome names from the EcoCyc'
        ' names to the bam file names. The names in the bam file should be '
        ' the same as the UCSC genome browser (they will be printed).')
    settings = parser.parse_args(argv)

    return settings

def get_coords(tfile):
    """
    Return the coords of genes in the file
    Arguments:
    - `tfile`: ptt.gz or rnt.gz file
    """
    coors = {}
    strand = {}
    tin = gzip.open(tfile)
    for line in tin:
        spl = line.strip().split()
        try:
            t1, t2 = spl[0].split('..')
            coors[spl[4]] = (int(t1)-1, int(t2))
            strand[spl[4]] = spl[1]
        except:
            pass
    return coors, strand

def main(argv=None):
    settings = process_command_line(argv)
    if len(settings.EC_chrlist) >= 2:
        chr_dict = dict(zip(
                settings.EC_chrlist.split(',')[0::2],
                settings.EC_chrlist.split(',')[1::2]))
    else:
        chr_dict = {}
    region_interactions, _, _, _=\
        RILseq.read_reads_table(open(settings.reads_in), settings.region)
    both_strs = defaultdict(lambda: defaultdict(int))
    if settings.summary:
        sig_reads = RILseq.read_significant_reads(
            settings.summary, chr_dict)
    if settings.known:
        known_reads = defaultdict(list)
        ptt_c, ptt_str = get_coords("%s.ptt.gz"%settings.refseq_dir)
        rnt_c, rnt_str = get_coords("%s.rnt.gz"%settings.refseq_dir)
        for line in open(settings.known):
            spl = line.strip().split()
            try:
                scoor = rnt_c[spl[0]]
                rcoor = ptt_c[spl[1]]
            except KeyError:
                pass
            else:
                for i in range(scoor[0], scoor[1]):
                    for j in range(rcoor[0], rcoor[1]):
                        known_reads[i].append(j)
                        known_reads[j].append(i)

    for reg1 in region_interactions:
        if reg1[2] != settings.chrn:
            continue
        for reg2 in region_interactions[reg1]:
            if reg2[2] != settings.chrn:
                continue
            if settings.summary:
                nsigs = 0
                for r1, r2, in region_interactions[reg1][reg2]:
                    nsigs += int((r2, reg2[1], reg2[2]) in \
                        sig_reads[(r1, reg1[1], reg1[2])])
            else:
                nsigs = len(region_interactions[reg1][reg2])
            both_strs[reg1[0]][reg2[0]] += nsigs
    if settings.sRNAs:
        from RILseq.ecocyc_parser import read_genes_data
        uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list, rRNAs = \
            read_genes_data(settings.ec_dir)
        sposs = set()
        for g in sRNAs_list:
            for i in range(uid_pos[g][1], uid_pos[g][2]):
                sposs.add(i)
    for r1 in both_strs:
        for r2 in both_strs[r1]:
            if both_strs[r1][r2] > settings.min_interactions:
                color = 'thickness=%dp'%max(
                    int(log(both_strs[r1][r2])/log(10)),1)
                if settings.sRNAs:
                    rset = set([i for i in range(r1, r1+settings.region)])
                    rset |= set([i for i in range(r2, r2+settings.region)])
                    if rset & sposs:
                        color = 'color=orange'
                if settings.known:
                    for k in set(range(r1, r1+settings.region)) & set(known_reads.keys()):
                        if set(range(r2, r2+settings.region)) & set(known_reads[k]):
                            color = 'color=red'
                            
                sys.stdout.write('%s %d %d %s %d %d %s\n'%(
                        settings.print_chr, r1+1, r1+settings.region,
                        settings.print_chr, r2+1, r2+settings.region,
                        color))

    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
