#!/usr/bin/env python

"""
Read the table of chimeric fragments and find regions that are over-represented
in their interactions. This script should remove non-specific chimeras using
Fisher excat's test by counting the number of fragments that support a pair of
genomic reginos and comparing to the expected number of interactions if this
interaction is random. The script treat the 5' end of the fragment and the
3' end separately since we observed that sRNAs tend to be in the 3' end
of the fragment.
Prior to testing, fragments that are probably not a consequence of ligation are
removed. This will be done by excluding reads that map less than 1000 bp apart
or on the same transcript (from EcoCyc annotations)
"""

import sys
import argparse
from collections import defaultdict

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
        description='Find over-represented regions of interactions.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'reads_in',
        help='An output file of map_chimeric_fragments.py with the chimeric'
        ' fragments.')
    parser.add_argument(
        '--ec_dir', #default='/home/users/assafp/Database/EcoCyc/19.0/data',
        help='EcoCyc data dir, used to map the regions to genes. If not'
        ' given only the regions will be reported')
    parser.add_argument(
        '--ribozero', default=False, action='store_true',
        help='Remove rRNA from the list of chimeric reads.')
    parser.add_argument(
        '--all_interactions', default=False, action='store_true',
        help='Skip all statistical tests and report all the interactions.')
    parser.add_argument(
        '--only_singles', default=False, action='store_true',
        help='Return only the single interactions. This should be used with'
        ' --all_interactions to count the number of single reads in the '
        'library.')
    parser.add_argument(
        '--est_utr_lens', type=int, default=100,
        help='Estimated UTRs lengths when there is not data.')
    parser.add_argument(
        '--EC_chrlist', default='chr,COLI-K12',
        help='A comma separated dictionary of chromosome names from the bam'
        ' file to EcoCyc names. See the names of chromosomes in bam file using'
        ' samtools view -H foo.bam.')
    parser.add_argument(
        '--refseq_dir', default='/home/users/assafp/EC/',
        help='RefSeq dir of organism to get the gene description from.')
    parser.add_argument(
        '-t', '--targets_file', #default='/home/users/assafp/Database/EcoCyc/19.0/data/curated_list_of_sRNA_mRNA_interactions.txt',
        help='A list of sRNA-mRNA interactions, should be in EcoCyc acc.')
    parser.add_argument(
        '-c', '--single_counts', #default='/home/hosts/disk19/E_coli/Sahar/FLAG_1-6_101-105_108_109_single_reads_counts_sum.txt',
        help='A file with the counts of single fragments per gene.')
    parser.add_argument(
        '-r', '--rep_table',
        help='A table containing data on REP elements. This file was generated'
        ' using SmartTables '
        '(e.g. this: http://ecocyc.org/group?id=biocyc14-8223-3640227683)')
    parser.add_argument(
        '-l', '--length', type=int, default=25,
        help='Length of sequence used for mapping. Used to determine the gene'
        ' in the regions.')
    parser.add_argument(
        '-s', '--shuffles', type=int, default=0,
        help='Shuffle the first sequence to compute an empirical p-value of the'
        ' hybridization energy using RNAup.')
    parser.add_argument(
        '--servers', #default=['pearl03']*16+['pearl01']*8+['gem']*8+['pearl02']*8+['jade']*8+['menash']*8+['amber']*10,
        help='A list of computers to run RNAup on (or number of CPUs')
    parser.add_argument(
        '--RNAup_cmd', default='RNAup -Xp -w 25 -b -o',
        help='Executable of RNAup with its parameters')
    parser.add_argument(
        '--pad_seqs', type=int, default=50,
        help='When computing RNAup pad the interacting regions.')
#    parser.add_argument(
#        '--maxdist', type=int, default=1000,
#        help='Maximal distance between mates to consider as same fragment.'
#        ' Here it will be used to reduce the counts in this distance from '
#        'either regions that interact with other regions of the genome.')
    parser.add_argument(
        '--seglen', type=int, default=100,
        help='Length of minimal segment of interaction.')
    parser.add_argument(
        '--maxsegs', type=int, default=5,
        help='Maximal number of consecutive segments that will be treated as'
        ' a region.')
    parser.add_argument(
        '--min_int', type=int, default=5,
        help='Minimal number of interactions to report.')
    parser.add_argument(
        '--max_pv', type=float, default=0.05,
        help='Maximal pvalue to report (after correction).')
    parser.add_argument(
        '--min_odds_ratio', type=float, default=2.0,
        help='Minimal odds ratio to report')
    



    settings = parser.parse_args(argv)

    return settings



def main(argv=None):
    settings = process_command_line(argv)
    if settings.ribozero:
        try:
            uid_pos,_,_,_,_,rRNAs = RILseq.ecocyc_parser.read_genes_data(
                settings.ec_dir)
        except IOError:
            rRNAs = None
        rr_pos = []
        chr_dict = dict(zip(
                settings.EC_chrlist.split(',')[1::2],
                settings.EC_chrlist.split(',')[0::2]))
        for rrgene in rRNAs:
            # Pad the position of the rRNA gene with the alignment lenght
            rr_pos.append([chr_dict[uid_pos[rrgene][0]]] +\
                              [uid_pos[rrgene][1]-settings.length] +\
                              [uid_pos[rrgene][2]+settings.length] +\
                              [uid_pos[rrgene][3]])
    else:
        rr_pos = None
    region_interactions, region_ints_as1, region_ints_as2, total_interactions=\
        RILseq.read_reads_table(
        open(settings.reads_in), settings.seglen, rr_pos, settings.only_singles)
    sys.stderr.write("Total interactions: %d\n"%total_interactions)

    # If all interactions are desired, skip the tests and report all
    if settings.all_interactions:
        interacting_regions = []
        for reg1, r1data in region_interactions.items():
            for reg2, clist in r1data.items():
                interacting_regions.append(
                    (1, len(clist), 0, reg1[0], reg1[0]+settings.seglen,
                     reg1[1], reg1[2], reg2[0], reg2[0]+settings.seglen,
                     reg2[1], reg2[2], 0, 0, 0))
    else:
    
        # Now run the test for each pair of interacting regions
        found_in_interaction = defaultdict(bool)
        interacting_regions = []
        # Start with the regions with the most interactions
        pairs_num = {}
        for reg1 in list(region_interactions.keys()):
            if region_ints_as1[reg1] < settings.min_int:
                continue
            for reg2 in list(region_interactions[reg1].keys()):
                if len(region_interactions[reg1][reg2]) >= settings.min_int:
                    pairs_num[(reg1, reg2)] = len(region_interactions[reg1][reg2])
        # Iterate the list of regions from the pairs with many interactions down
        for (reg1, reg2) in sorted(pairs_num, key=pairs_num.get, reverse=True):
            pv, ints, odds, r1_from, r1_to, r2_from, r2_to, mat_b, mat_c,mat_d=\
                RILseq.minpv_regions(
                reg1, reg2, region_interactions, region_ints_as1,
                region_ints_as2, total_interactions, found_in_interaction,
                settings.seglen, settings.maxsegs, settings.min_odds_ratio)
            pv *= len(pairs_num)
            if pv <= settings.max_pv:
                # Mark as participating
                for r1 in range(r1_from, r1_to, settings.seglen):
                    for r2 in range(r2_from, r2_to, settings.seglen):
                        found_in_interaction[
                            (r1, reg1[1], reg1[2], r2, reg2[1], reg2[2])] = True
                # Report this interaction
                interacting_regions.append(
                    (pv, ints, odds, r1_from, r1_to, reg1[1], reg1[2],  r2_from,
                     r2_to, reg2[1], reg2[2],  mat_b, mat_c, mat_d))
    # Read the additional data to decorate the results with
    RILseq.report_interactions(
        region_interactions, sys.stdout, interacting_regions, settings.seglen,
        settings.ec_dir, settings.EC_chrlist, settings.refseq_dir,
        settings.targets_file, settings.rep_table, settings.single_counts,
        settings.shuffles, settings.RNAup_cmd, settings.servers,
        settings.length, settings.est_utr_lens, settings.pad_seqs)
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
