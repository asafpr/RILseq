#!/usr/bin/env python

"""
Generate a gff file of transcripts using an EcoCyc data flat-files
"""

import sys
import argparse

from RILseq import ecocyc_parser

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(description='Generate EcoCyc genes gff file')

    parser.add_argument(
        'ec_dir',
        help='EcoCyc flat-files directory.')
    parser.add_argument(
        '--EC_chrlist', default='COLI-K12,chr',
        help='A comma separated dictionary of chromosome names from the '
        'EcoCyc name to the bam name. See the names of chromosomes in bam file'
        ' using samtools view -H foo.bam.')
    parser.add_argument(
        '--est_utr_lens', type=int, default=100,
        help='Estimated UTRs lengths when there is not data.')
    settings = parser.parse_args(argv)

    return settings

def main(argv=None):
    settings = process_command_line(argv)
    if settings.EC_chrlist:
        chr_dict = dict(zip(
                settings.EC_chrlist.split(',')[0::2],
                settings.EC_chrlist.split(',')[1::2]))
    else:
        chr_dict = None
    ecocyc_parser.generate_transcripts_file(
        sys.stdout, utr_len=settings.est_utr_lens,
        ec_dir=settings.ec_dir, chr_dict=chr_dict)
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
