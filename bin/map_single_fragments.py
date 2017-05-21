#!/usr/bin/env python

"""
This script is used to map fastq files to the genome. The input is a comma
separated list of fastq[.gz] files (or two lists if the input is paired-end).
The output are bam files with the mapped reads, a table containing the number
of reads mapped to each gene and a wiggle file with the coverage of the
reads.
"""

import sys
import argparse
import pysam
import csv
import os

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
        description='Map fastq files to the genome using bwa.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genome_fasta',
        help='Name of genome fasta file. The file must be indexed using'
        'bwa index command prior to this run.')
    parser.add_argument(
        '-1', '--fastq_1', action='append', nargs='+',
        help='A list of the first read of the sequencing.')
    parser.add_argument(
        '-2', '--fastq_2', action='append', nargs='*',
        help='A list of the second read of the sequencing.'
        ' The order of these files should be as same as -1. Optional.')
    parser.add_argument(
        '-g', '--genes_gff',
        help='Name of gff file to count the reads per gene. If not given '
        ' or not readable, skip this stage.')
    parser.add_argument(
        '-r', '--reverse_complement', default=False,
        action='store_true',
        help='Treat the reads as reverse complement only when counting'
        ' number of reads per gene and generating wig file. The resulting BAM'
        ' files will be the original ones. Use this when treating libraries'
        " built using Livny's protocol.")
    parser.add_argument(
        '-f', '--feature', default='exon',
        help='Name of features to count on the GTF file (column 2).')
    parser.add_argument(
        '-i', '--identifier', default='gene_id',
        help='Name of identifier to print (in column 8 of the GTF file).')
    parser.add_argument(
        '-v', '--overlap', type=int, default=5,
        help='Minimal required overlap between the fragment and the feature.')
    parser.add_argument(
        '-m', '--allowed_mismatches', type=float, default=2,
        help="Allowed mismatches for BWA mapping.")
    parser.add_argument(
        '-o', '--outhead', default='bwa_mapped_single_reads',
        help='Output file names of counts table (suffixed _counts.txt) and'
        ' wiggle file (suffixed _coverage.wig)')
    parser.add_argument(
        '-d', '--dirout', default='.',
        help='Output directory, default is this directory.')
    parser.add_argument(
        '--bwa_exec', default='bwa',
        help='bwa command')
    parser.add_argument(
        '-S', '--samtools_cmd', default='samtools',
        help='Samtools executable.')
    parser.add_argument(
        '-a', '--params_aln', default='-t 32 -R 200',
        help='Additional parameters for aln function of bwa.')
    parser.add_argument(
        '-s', '--sampe_params', default='-a 1500 -P',
        help='Additional parameters for sampe function of bwa.')
    parser.add_argument(
        '--samse_params', default=' ',
        help='Additional parameters for samse function of bwa.')

    settings = parser.parse_args(argv)

    return settings

def main(argv=None):
    settings = process_command_line(argv)
    if not os.path.exists(settings.dirout):
        os.makedirs(settings.dirout)
    outwig = open("%s/%s_coverage.wig"%(settings.dirout, settings.outhead), 'w')
    if settings.genes_gff:
        try:
            pos_feat_list, all_features = RILseq.read_gtf(
                open(settings.genes_gff), settings.feature, settings.identifier)
        except IOError:
            settings.genes_gff = None
        gcounts = {}
        lib_order = []
    fastq_2_list = list(RILseq.flat_list(settings.fastq_2))
    for i, r1_name in enumerate(RILseq.flat_list(settings.fastq_1)):
        try:
            r2_name = fastq_2_list[i]
        except IndexError:
            r2_name = None
        outhead = r1_name.rsplit('.', 1)[0]
        libname = outhead.rsplit('/',1)[-1]
        outhead = '%s_bwa'%libname
        bamname = RILseq.run_bwa(
            settings.bwa_exec, r1_name, r2_name,
            settings.dirout, outhead, settings.allowed_mismatches,
            settings.genome_fasta, settings.params_aln, settings.sampe_params,
            settings.samse_params, settings.samtools_cmd)
        samfile = pysam.Samfile(bamname)
        if settings.genes_gff:
            lib_order.append(libname)
            gcounts[libname] = RILseq.count_features(
                pos_feat_list, samfile, settings.overlap,
                rev=settings.reverse_complement)
        coverage = RILseq.generate_wig(
            samfile, rev=settings.reverse_complement, first_pos=False)
        RILseq.print_wiggle(
            coverage, "%s_single_fragments_coverage"%libname,
            "%s single fragments coverage"%libname, outwig)
    # Print the table of counts
    if settings.genes_gff:
        outtable = open(
            "%s/%s_counts.txt"%(settings.dirout, settings.outhead), 'w')
        outt = csv.writer(outtable, delimiter='\t')
        outt.writerow(['Gene name'] + lib_order)
        for g in sorted(list(all_features)):
            row_out = [g]
            for libn in lib_order:
                row_out.append(gcounts[libn][g])
            outt.writerow(row_out)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
