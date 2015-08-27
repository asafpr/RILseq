#!/usr/bin/env python

"""
Given already mapped fusions using the reads file (format:
gene1 gene2 position1 strand1 position2 strand2 read_name)
Use the original BAM to plot the ends of the reads as BED file to be presented
by the genome browser.
Color code as specified in the parametrs
"""

import sys
import argparse
import csv
from collections import defaultdict

from Bio.Seq import Seq
import pysam
from Bio import SeqIO

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(
        description='Generate BED graph of the reads.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'genome',
        help='genome fasta file.')
    parser.add_argument(
        'list_reads',
        help='File with list of reads and their fused positions.')
    parser.add_argument(
        'bamfile',
        help='The original bam file with the full reads.')
    parser.add_argument(
        'track_name',
        help='Name of track')
    parser.add_argument(
        'track_desc',
        help='Description of the track')


    parser.add_argument(
        '-r', '--reverse', default=False, action='store_true',
        help='The original bam file is the reverse complement of the RNA.')
    parser.add_argument(
        '-s', '--summary',
        help='Print only reads that are found to be significant in this summary file.')
    parser.add_argument(
        '-e', '--gene_name',
        help='Print reads involve only this gene (EcoCyc ID), '
        'applies only with -s')
    parser.add_argument(
        '--pos_first', default='255,0,0',
        help='Color of first part, positive strand.')
    parser.add_argument(
        '--pos_second', default='51,102,255',
        help='Color of second part, positive strand.')
    parser.add_argument(
        '--rev_first', default='255,0,0',
        help='Color of first part, reverse strand.')
    parser.add_argument(
        '--rev_second', default='51,102,255',
        help='Color of second part, reverse strand.')
    parser.add_argument(
        '--EC_chrlist', default='COLI-K12,chr',
        help='A comma separated dictionary of chromosome names from the EcoCyc'
        ' names to the bam file names. The names in the bam file should be '
        ' the same as the UCSC genome browser (they will be printed).')
    settings = parser.parse_args(argv)
    return settings

def get_reads_seqs(bamfile, rnames, rev=False):
    """
    Return the sequences of all the reads from the bam file
    Arguments:
    - `bamfile`: The pysam file
    - `rnames`: reads names
    """
    r1_seqs = {}
    r2_seqs = {}
    rqns = set()
    reads = defaultdict(list)
    for read in bamfile.fetch(until_eof=True):
        rqns.add(read.qname)
        reads[read.qname].append(read)
    for rn in set(rnames) & rqns:
        for read in reads[rn]:
            if read.is_read1==rev:
                outseq = Seq(read.seq)
                if not read.is_reverse:
                    outseq = outseq.reverse_complement()
                r1_seqs[read.qname] = str(outseq)
            else:
                outseq = Seq(read.seq)
                if read.is_reverse:
                    outseq = outseq.reverse_complement()
                r2_seqs[read.qname] = str(outseq)
    # r1_seqs is the 3' end of the second fused RNA, r2_seqs is the 5' of the
    # first fused RNA
    return r1_seqs, r2_seqs

def extend_alignment(rseq, pos5p, pos3p, is_read1, strand, genome, mismatch=1):
    """
    Align the rseq to the genome in the specified position. Return the last
    position of the read mapped to the genome.
    Use local alignment
    Arguments:
    - `rseq`: Read sequence
    - `pos5p`: the 5' position, exact if read 2 or as limit if read 1
    - `pos3p`: the 3' position, exact if read 1 or as limit if read 2
    - `is_read1`: This read is read 1
    - `strand`: mapping strand
    - `genome`: The genome Seq object
    - `mismatch`: allowed mismatches
    """
    rcnt = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    glen = len(genome)
    if is_read1:
        # Start from the last position and move on to the 5' end
        if strand == '-':
            ipos = 0
            for _ in range(mismatch+1):
                try:
                    while rcnt[genome[(pos3p+ipos)%glen]] == rseq[-(ipos+1)]:
                        ipos += 1
                except IndexError:
                    return ipos - 1
                ipos += 1
            return ipos
        else:
            ipos = 0
            for _ in range(mismatch+1):
                try:
                    while genome[(pos3p-ipos)%glen] == rseq[-(ipos+1)]:
                        ipos += 1
                except IndexError:
                    return ipos-1
                ipos += 1
            return ipos
    else:
        if strand == '-':
            ipos = 0
            for _ in range(mismatch+1):
                try:
                    while rcnt[genome[(pos5p-ipos)%glen]] == rseq[ipos]:
                        ipos += 1
                except IndexError:
                    return ipos -1
                ipos += 1
            return ipos 
        else:
            ipos = 0
            for _ in range(mismatch+1):
                try:
                    while genome[(pos5p+ipos)%glen] == rseq[ipos]:
                        ipos += 1
                except IndexError:
                    return ipos - 1
                ipos += 1
            return ipos 
        
                
def read_significant_reads(summary_file, chr_dict, gname=None):
    """
    Return a dict from r1->[r2] of regions that are significant
    Arguments:
    - `summary_file`: A summary results file
    - `chr_dict`: Dictionary from chr name in summary file (EcoCyc) to bam
    - `gname`: Choose reads of only this gene
    """
    sig_reg = defaultdict(list)
    for line in csv.DictReader(open(summary_file), delimiter='\t'):
        r1_from = int(line['RNA1 from'])-1
        r1_to = int(line['RNA1 to'])
        try:
            r1_chrn = chr_dict[line['RNA1 chromosome']]
        except KeyError:
            r1_chrn = line['RNA1 chromosome']
        r1_str = line['RNA1 strand']
        r2_from = int(line['RNA2 from'])-1
        r2_to = int(line['RNA2 to'])
        try:
            r2_chrn = chr_dict[line['RNA2 chromosome']]
        except KeyError:
            r2_chrn = line['RNA2 chromosome']
        r2_str = line['RNA2 strand']
        
        if gname:
            try:
                if (gname not in line['RNA1 EcoCyc ID']) and\
                        (gname not in line['RNA2 EcoCyc ID']):
                    continue
            except KeyError:
                continue
        for i in range(r1_from, r1_to):
            for j in range(r2_from, r2_to):
                sig_reg[(i, r1_str, r1_chrn)].append((j, r2_str, r2_chrn))
    return sig_reg
        
def find_overlap(s1, s2):
    """
    Find overlaps between two reads. Assume they are both in the same
    orientation (r1 is revcomp)
    Return 3 seuqnces: s1, overlap, s2
    Arguments:
    - `s1`: first sequence, this is mate 2 actually in our experiments
    - `s2`: last sequence, mate 1
    """
    for i in range(min(len(s1), len(s2)))[::-1]:
        if s1[-i:]==s2[:i]:
            return s1[:-i], s1[-i:], s2[i:]
    return s1, '', s2


def main(argv=None):
    settings = process_command_line(argv)
    # Read the read names and positions
    read_5ps = {}
    read_3ps = {}
    read_genes = {}
    genome = {}
    gsize = {}
    for sr in SeqIO.parse(settings.genome, 'fasta'):
        genome[sr.id] = sr.seq
        gsize[sr.id] = len(sr.seq)
    if len(settings.EC_chrlist) >= 2:
        chr_dict = dict(zip(
                settings.EC_chrlist.split(',')[0::2],
                settings.EC_chrlist.split(',')[1::2]))
    else:
        chr_dict = {}
    if settings.summary:
        sig_reads = read_significant_reads(
            settings.summary, chr_dict, gname=settings.gene_name)

    for line in csv.reader(open(settings.list_reads), delimiter='\t'):
        # skip single
        if len(line) > 7 and line[7]=="single":
            continue
        if settings.summary:
            if (int(line[4])-1, line[5], line[3]) not in\
                    sig_reads[(int(line[1])-1, line[2], line[0])]:
                continue
        read_5ps[line[6]] = [int(line[1])-1, line[2], line[0]]
        read_3ps[line[6]] = [int(line[4])-1, line[5], line[3]]
#        read_genes[line[6]] = [line[0], line[1]]
    # Read the bam file and return the long sequences
    r1_seqs, r2_seqs = get_reads_seqs(
        pysam.Samfile(settings.bamfile), read_5ps.keys(), rev=settings.reverse)
    # For each read find the overlap, if exists and find the fusion point
    outer = csv.writer(sys.stdout, delimiter='\t')
    print 'track name="%s" description="%s" visibility=4 itemRgb="On"'%(
        settings.track_name, settings.track_desc)
    # Because I'm lazy, the code is written so r1 is the 3' end of the fragment
    for rname in set(r2_seqs.keys()):
        if rname in r1_seqs:
            r2seq = r2_seqs[rname]
            r1seq = r1_seqs[rname]
        else: # single-end
            r2seq = r2_seqs[rname]
            r1seq = ''
        s1, overlap, s2 = find_overlap(r2seq, r1seq)
        side_5p_len = extend_alignment(
            s1+overlap+s2, read_5ps[rname][0], 0, False, read_5ps[rname][1],
            genome[read_5ps[rname][2]])
        side_3p_len = extend_alignment(
            s1+overlap+s2, 0, read_3ps[rname][0], True, read_3ps[rname][1],
            genome[read_3ps[rname][2]])
        # Write each of the sides to the output file
        if read_5ps[rname][1] == '+':
            gfrom = max(0,  read_5ps[rname][0])
            gto = min(gsize[read_5ps[rname][2]], read_5ps[rname][0]+side_5p_len)
            outer.writerow([
                    read_5ps[rname][2], gfrom, gto, "%s_5p"%rname, 0, '+',
                    gfrom, gto, settings.pos_first])
        elif read_5ps[rname][1] == '-':
            gfrom = max(0, read_5ps[rname][0]-side_5p_len+1)
            gto = min(gsize[read_5ps[rname][2]], read_5ps[rname][0]+1)
            outer.writerow([
                    read_5ps[rname][2], gfrom, gto, "%s_5p"%rname, 0, '-',
                    gfrom, gto,settings.rev_first])
        if read_3ps[rname][1] == '+':
            gfrom = max(0, read_3ps[rname][0]-side_3p_len+1)
            gto = min(gsize[read_3ps[rname][2]], read_3ps[rname][0]+1)
            outer.writerow([
                    read_3ps[rname][2], gfrom, gto,"%s_3p"%rname, 0, '+',
                    gfrom, gto, settings.pos_second])
        elif read_3ps[rname][1] == '-':
            gfrom = max(0, read_3ps[rname][0])
            gto = min(gsize[read_3ps[rname][2]], read_3ps[rname][0]+side_3p_len)
            outer.writerow([
                    read_3ps[rname][2], gfrom, gto, "%s_3p"%rname, 0, '-',
                    gfrom, gto, settings.rev_second])
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
