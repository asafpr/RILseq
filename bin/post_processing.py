"""
This script is used to run post processing analysis on the RILseq results.
This will generate table S1 in the RILseq essay.
"""

import argparse
import pysam
import csv
import glob
import os
import sys
import gzip
from subprocess import check_output, Popen, PIPE


def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    :param argv: is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Post processing analysis of RILseq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'process_nextseq_work_dir',
        help='Name of the process_nextseq output directory containing fastq files and BAM files')
    parser.add_argument(
        '-R', '--RILseq_work_dir', default=None,
        help='Name of the RILseq work directory containing RILseq output files.')
    parser.add_argument(
        '--count_path', default='/cs/icore/asafp5/lab_nextseq_complete_data/bin/ProcessNextSeq/count_PE_fragments.py',
        help='path for the count_PE_fragments script for counting')

    parser.add_argument(
        '-g', '--gtf',
        help='GTF file containing the features.')
    parser.add_argument(
        '-r', '--reverse', action='store_true', default=False,
        help='Count the genes on the reverse strand of the mapping.')
    parser.add_argument(
        '-f', '--feature', default='exon',
        help='Name of features to count on the GTF file (column 2).')
    parser.add_argument(
        '-i', '--identifier', default='gene_id',
        help='Name of identifier to print (in column 8 of the GTF file).')
    parser.add_argument(
        '-o', '--overlap', type=int, default=5,
        help='Minimal required overlap between the fragment and the feature.')
    parser.add_argument(
        '-c', '--skip_clipped', action='store_true', default=False,
        help='Discard reads with soft clipped bases, when generating counts file.')
    parser.add_argument(
        '--output_table', default='S1.txt',
        help='Output table for writing of the S1 table.')

    settings = parser.parse_args(argv)
    return settings


def generate_table_s1(settings, outfile):
    """
    the function generates the S1 table using the work directory given and other settings.
    :param settings: se
    :param outfile:
    :return:
    """
    outer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
    header = ['Library name', 'Total # of reads', '# of reads after processing', 'Processed reads/Total reads',
              '# of mapped reads (strict)', 'Multiple read counts', 'Antisense read counts', 'IGR read counts',
              '# of mapped reads (non strict)', 'Mapped reads/Processed reads(strict)',
              'Mapped reads/Processed reads(non strict)']
    if settings.RILseq_work_dir:
        header.extend(['Total number of fragments', 'Mapped as single', 'Mapped as single excluding rRNAs',
                       'Mapped as chimeric', 'Mapped as chimeric excluding rRNAs', '% of Mapped fragments',
                       '% of Chimeric fragments out of mapped fragments',
                       '# of statistically significant interactions (S-chimeras)',
                       '# of fragments in S-chimera', '% of fragments in S-chimeras out of chimeric fragments',
                       '# of mapped reads (original)'])

    content = []
    no_fastq = False

    if len(glob.glob("%s/fastq_files/*_1.fastq" % settings.process_nextseq_work_dir)) == 0:
        if len(glob.glob("%s/fastq_files/*_1.fastq.gz" % settings.process_nextseq_work_dir)) == 0:
            print ("There are no fastq files in the given process_next_seq directory, Exiting")
            return
        else:
            file_list = glob.glob("%s/fastq_files/*_1.fastq.gz" % settings.process_nextseq_work_dir)
            no_fastq = True

    else:
        file_list = glob.glob("%s/fastq_files/*_1.fastq" % settings.process_nextseq_work_dir)

    for fastqfile in file_list:
        if not no_fastq:
            libname = fastqfile.split('_1.fastq')[0]
            total_reads = sum(1 for line in open(fastqfile) if line.startswith('@'))
            if os.path.exists(fastqfile.replace('_1.fastq', '_2.fastq')):
                total_reads *= 2

        else:  # there are only fastq.gz files and no fastq files
            libname = fastqfile.split('_cutadapt_1.fastq.gz')[0]
            total_reads = 0

        procfile = libname + '_cutadapt_1.fastq.gz'
        processed_reads = sum(1 for line in gzip.open(procfile) if line.startswith('@'))
        if os.path.exists(procfile.replace('_1.fastq.gz', '_2.fastq.gz')):
            processed_reads *= 2
        if total_reads > 0:
            processed_fraction = "{:1.2f}".format(processed_reads/float(total_reads))
        else:
            processed_fraction = "0"
        bamfile = libname.replace('fastq_files', 'bwa_mapping')+'_cutadapt_bwa.bam'

        count_cmd = ['python', settings.count_path, '-g', settings.gtf, '-s', bamfile, '-o', str(settings.overlap),
                     '-f', settings.feature, '-i', settings.identifier, '--sum']
        if settings.reverse:
            count_cmd.append('-r')
        if settings.skip_clipped:
            count_cmd.append('-c')
        count_cmd_strict = count_cmd + ['-t']

        print (count_cmd)
        mapped_reads = Popen(' '.join(count_cmd), shell=True, stdout=PIPE).stdout.read().strip()
        mapped_reads_strict, multiple, antisense, igr = Popen(' '.join(count_cmd_strict), shell=True, stdout=PIPE)\
            .stdout.read().strip().split(',')

        mapped_fraction = "{:1.2f}".format(int(mapped_reads)/float(processed_reads))
        mapped_fraction_strict = "{:1.2f}".format(int(mapped_reads_strict)/float(processed_reads))

        row = [libname.split('/')[-1], str(total_reads), str(processed_reads), processed_fraction,
                        mapped_reads_strict, multiple, antisense, igr, mapped_reads, mapped_fraction_strict,
                        mapped_fraction]

        if settings.RILseq_work_dir:
            short_name = libname.split('/')[-1]
            fastq1 = settings.RILseq_work_dir+short_name+"_cutadapt_bwa_ends_1.fastq"
            fastq2 = settings.RILseq_work_dir+short_name+"_cutadapt_bwa_ends_2.fastq"
            if os.path.exists(fastq1):
                numfrag = len(list(set([line.strip() for line in open(fastq1) if line.startswith('@')]).intersection(
                        [line.strip() for line in open(fastq2) if line.startswith('@')])))

                frags_file = settings.RILseq_work_dir+short_name+'_cutadapt_bwa.bam_all_fragments_l25.txt'
                singles = sum(1 for line in open(frags_file) if "single" in line)
                chimeras = sum(1 for line in open(frags_file) if "chimera" in line)

                singles_file = settings.RILseq_work_dir+short_name+'_cutadapt_bwa.bam_all_fragments_l25.txt_only_single.txt'
                all_file = settings.RILseq_work_dir+short_name+'_cutadapt_bwa.bam_all_fragments_l25.txt_all_interactions.txt'
                sig_file = settings.RILseq_work_dir+short_name+'_cutadapt_bwa.bam_all_fragments_l25.txt_sig_interactions_with-total.txt'

                temp = open(singles_file)
                temp.next()
                map_si_exc_rrna = sum(int(line.split('\t')[14]) for line in temp)

                temp = open(all_file)
                temp.next()
                map_ch_exc_rrna = sum(int(line.split('\t')[14]) for line in temp)

                per_mapped = "{:1.1f}%".format(((singles+chimeras)/float(numfrag))*100)
                per_chim = "{:1.1f}%".format((chimeras/float(chimeras+singles))*100)

                stat_interactions = sum(1 for line in open(sig_file)) - 1
                temp = open(sig_file)
                temp.next()
                stat_frags = sum(int(line.split('\t')[14]) for line in temp)
                per_stat = "{:1.1f}%".format((stat_frags/float(chimeras))*100)

                row.extend([numfrag, singles, map_si_exc_rrna, chimeras, map_ch_exc_rrna, str(per_mapped), str(per_chim),
                            stat_interactions, stat_frags, per_stat])
        countsfile = libname.replace('fastq_files', 'bwa_mapping')+'_cutadapt_bwa.counts'
        if os.path.exists(countsfile):
            counts = sum(int(line.split('\t')[1]) for line in open(countsfile))
            row.extend([str(counts)])

        content.append(row)

    outer.writerow(header)
    for lib in content:
        outer.writerow(lib)


def check_read(read):
    """
    a callback function for pysam.count function, returns true if the read should be included in the count
    :param read: the read to check
    :return: True of False whether the read should be included in the count
    """
    return read.is_proper_pair  # and not (read.is_secondary or read.is_supplementary) \
#        and read.mapping_quality >= 10


def main(argv=None):
    settings = process_command_line(argv)
    generate_table_s1(settings, open(settings.output_table, 'w'))

    return 0  # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
