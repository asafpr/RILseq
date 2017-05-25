from subprocess import call
import logging
from tempfile import NamedTemporaryFile
from collections import defaultdict
import csv
import sys
from Bio.Seq import Seq
import numpy as np
from scipy.stats import fisher_exact
import glob
import gzip


from . import ecocyc_parser
from . import RNAup_tar_pred
from . import dinuc_shuffle

# A small helper function to handle argmunets parsing
def flat_list(list_to_flat):
    """
    Flattent the list of arguments to one list. The lists of input might be list
    of lists, this function yields them as one list
    Arguments:
    - `list_to_flat`: A nested list
    """
    if not isinstance(list_to_flat, list):
        yield list_to_flat
    else:
        for item in list_to_flat:
            for it2 in flat_list(item):
                yield(it2)

# Functions for simple mapping (single fragment)
def run_bwa(bwa_cmd, fname1, fname2, output_dir, output_prefix, mismatches,
            fasta_genome, params_aln, params_sampe, params_samse, samtools_cmd):
    """
    Run bwa on paired or single end fastq files. write a sorted bam and a bai
    file
    Arguments:
    - `bwa_cmd`: executable of bwa
    - `fname1`: The R1 fastq file
    - `fname2`: The R2 fastq file
    - `output_dir`: Where to write the files to
    - `output_prefix`: Name of bam file (without the extension
    - `mismatches`: Allowed mismatches
    - `fasta_genome`: genome fasta file. Must be indexed!
    - `params_aln`: extra parametrs for aln command
    - `params_sampe`: extra paarmeters for sampe execution
    - `params_samse`: extra parameters for samse execution
    - `samtools_cmd`: samtools command

    Return
    - `bamfile`: Output bam file name
    """
    # Run aln of bwa on both files
    sai1 = NamedTemporaryFile(dir=output_dir)
   
    # Added 11.2.17 - bug fix for adding of the params_aln parameters.
    # This param was ignored when using the splitting option of subprocess.
    next_cmd = [bwa_cmd, 'aln', '-n', str(mismatches)]
    next_cmd.extend(params_aln.split())
    next_cmd.extend([fasta_genome, fname1])

    logging.info("Executing %s", ' '.join(next_cmd))
    call(next_cmd, stdout=sai1)

    if fname2:
        sai2 = NamedTemporaryFile(dir=output_dir)
        next_cmd = [bwa_cmd, 'aln', '-n', str(mismatches)]
        next_cmd.extend(params_aln.split())
        next_cmd.extend([fasta_genome, fname2])
        logging.info("Executing %s", ' '.join(next_cmd))
        call(next_cmd, stdout=sai2)
        # Run sampe on both sai files

    # Changed 11.2.17 by niv. Piping bug fix in version update of the cluster
    bwa_file = NamedTemporaryFile(dir=output_dir)
    if fname2:
        bwa_sam_cmd = ' '.join([bwa_cmd, 'sampe', params_sampe, fasta_genome, sai1.name, sai2.name, fname1, fname2,
                                '>', bwa_file.name])
    else:  # Single end
        bwa_sam_cmd = ' '.join([bwa_cmd, 'samse', params_samse, fasta_genome, sai1.name, fname1, '>', bwa_file.name])

    view_file = NamedTemporaryFile(dir=output_dir)
    view_cmd = ' '.join([samtools_cmd, 'view', '-u', bwa_file.name, '>', view_file.name])
    sort_cmd = ' '.join([samtools_cmd, 'sort', view_file.name, '-o', "%s/%s.bam" % (output_dir, output_prefix)])

    logging.info("Executing %s" % bwa_sam_cmd)
    call(bwa_sam_cmd, shell=True)
    logging.info("Executing %s" % view_cmd)
    call(view_cmd, shell=True)
    logging.info("Executing %s" % sort_cmd)
    call(sort_cmd, shell=True)

    bamname = "%s/%s.bam"%(output_dir, output_prefix)
    # Indexing the bam file
    index_cmd = [samtools_cmd, 'index', bamname]
    logging.info("Indexing bam file %s"%' '.join(index_cmd))
    call(index_cmd)
    return bamname



def read_gtf(gtf_file, feature, identifier):
    """
    Read a GTF file and return a list in the length of the genome in which each
    position contains a list of features that overlap this position
    Arguments:
    - `gtf_file`: An open gtf_file
    - `feature`: the name of the feature to index
    - `identifier`: The identifier to use from column 8
    """
    # First initialize a dictionary
    # posfeat->[chromosome+strand]->[position]->
    # [set of entities in this position]
    pos_feat = defaultdict(lambda: defaultdict(set))
    # Get all the names of the features
    all_features = set(['~~intergenic', '~~antisense'])
    for line in csv.reader(gtf_file, delimiter='\t'):
        if line[2] != feature:
            continue
        ids_dict = {}
        for id_pair in line[8].strip().split(';'):
            try:
                k, v = id_pair.strip().split(' ')
            except ValueError:
                pass
            ids_dict[k] = v.replace('"','')
        fid = ids_dict[identifier]
        all_features.add(fid)
        # Change to 0-based coordinates and add this feature to all the
        # positions it convers
        for i in range(int(line[3])-1, int(line[4])):
            # Concat the strand to the name of the chromosome
            pos_feat[line[0]+line[6]][i].add(fid)
    # Change the dictionary to list
    pos_feat_list = {}
    for chrom, data in pos_feat.items():
        maxpos = max(data.keys())
        list_of_sets = []
        for k in range(maxpos+1):
            list_of_sets.append(list(data[k]))
        pos_feat_list[chrom] = list_of_sets
    return pos_feat_list, all_features



def get_paired_pos(read, rev=False):
    """
    Given a read which is the positive of the pairs return a strand and
    start and end positions
    Arguments:
    - `read`: A read from pysam
    - `rev`: The read is the reverse complement of the RNA, in this case
             return the opposite strand
    """
    strand = '+'
    if rev!=read.is_read2:
        strand = '-'
    fpos = read.pos
    tpos = read.tlen + fpos
    return strand, fpos, tpos

def get_single_pos(read, rev=False):
    """
    Given a read which is the positive of the pairs return a strand and
    start and end positions
    Arguments:
    - `read`: A read from pysam
    - `rev`: The read is the reverse complement of the RNA, in this case
             return the opposite strand
    """
    strand = '+'
    if rev!=read.is_reverse:
        strand = '-'
    fpos = read.pos
    tpos = read.qlen + fpos
    return strand, fpos, tpos




def count_features(
    features_lists, samfile, overlap, rev=False, checkpoint=1000000,
    get_sum=False):
    """
    Go over the samfile and for each pair of reads find the features that
    overlap the fragment with at least 'overlap' nucleotides. Add 1 to the count
    of these features
    
    Arguments:
    - `features_lists`: The list of features returned from the read_gtf function
    - `samfile`: A pysam object
    - `overlap`: The minimal overlap between the feature and read
    - `rev`: reverse the strand of the read
    - `checkpoint`: Report every 100000 reads processed, set to None or False
                    for silencing
    - `get_sum`: Return the number of reads as well
                    
    Return:
    - `fcounts`: A dictionary from gene name to number of reads mapped to the
                 gene. ~~antisense and ~~intergenic count the number of reads
                 that were not mapped to a gene and found antisense to a gene
                 or in intergenic region.
    - `reasd_num`: Number of reads if get_sum is True
    """
    fcounts = defaultdict(int)
    counter = 0
    for read in samfile.fetch():
        if read.is_paired:
            if read.is_reverse or read.is_unmapped or\
                read.mate_is_unmapped or\
                read.is_reverse==read.mate_is_reverse or\
                not read.is_proper_pair:

                continue
        else: # single end
            if  read.is_unmapped:
                continue
        # Take only the forward mate
        counter += 1
        if checkpoint and counter%checkpoint==0:
#            pass
            sys.stderr.write("Processed %i fragments\n"%counter)
        try:
            chrname = samfile.getrname(read.tid)
        except ValueError:
            sys.stderr.write(str(read)+"\n")
        # Get the positions of the fragment
        if read.is_paired:
            strand, fpos, tpos = get_paired_pos(read, rev=rev)
        else:
            strand, fpos, tpos = get_single_pos(read, rev=rev)
        
        # Count the number of times a feature intersects with the fragment
        rcounts = defaultdict(int)
        for fset in features_lists[chrname+strand][fpos:tpos]:
            for el in fset:
                rcounts[el] += 1
        # Go over the list of features, if the number of counts is above the
        # Threshold add 1 to the count of this feature
        is_counted = False
        for feature, counts in rcounts.items():
            if counts >= overlap:
                fcounts[feature] += 1
                is_counted = True
        if not is_counted:
            # Test if antisense
            rev_str = '-'
            if strand == '-':
                rev_str = '+'
            rev_counts = defaultdict(int)
            for fset in features_lists[chrname+rev_str][fpos:tpos]:
                for el in fset:
                    rev_counts[el] += 1
            is_antis = False
            for feature, counts in rev_counts.items():
                if counts >= overlap:
                    is_antis = True
                    break
            if is_antis:
                fcounts['~~antisense'] += 1
            else:
                fcounts['~~intergenic'] += 1
    if get_sum:
        return fcounts, counter
    else:
        return fcounts


def generate_wig(samfile, rev=False, first_pos=False):
    """
    Go over the samfile and return two histograms (for + and - strands) of
    coverage
    
    Arguments:
    - `samfile`: A pysam object
    - `rev`: reverse the strand of the read
    - `first_pos`: Count only the first position of each read
    """
    # Build the structure of the dictionary chromosome->strand->list of 0
    coverage = {}
    for i, rfg in enumerate(samfile.references):
        rlen = samfile.lengths[i]
        coverage[rfg] = {'-':[0] * rlen, '+':[0] * rlen}
    for read in samfile.fetch():
        if read.is_paired:
            if read.is_reverse or read.is_unmapped or\
                read.mate_is_unmapped or\
                read.is_reverse==read.mate_is_reverse or\
                not read.is_proper_pair:
                continue
        else: # single end
            if  read.is_unmapped:
                continue
        # Take only the forward mate
        try:
            chrname = samfile.getrname(read.tid)
        except ValueError:
            logging.warn("Read has no valid chr name %s"%(str(read)))
            continue
        # Get the positions of the fragment
        if read.is_paired:
            strand, fpos, tpos = get_paired_pos(read, rev=rev)
        else:
            strand, fpos, tpos = get_single_pos(read, rev=rev)
        rrange = range(fpos, tpos)
        if first_pos:
            if strand == '+':
                rrange = [fpos]
            else:
                rrange = [tpos]
        for i in rrange:
            try:
                coverage[chrname][strand][i] += 1
            except IndexError:
                logging.warn("IndexError: trying to set index %d on chr %s, bu length is only %d"%(i, chrname, len(coverage[chrname][strand])))
    return coverage


def print_wiggle(coverage, title, description, outf):
    """
    Print the coverage into an open wiggle file
    Arguments:
    - `coverage`: returned from generate_wig
    - `title`: Title of wig track
    - `description`: Description of track
    - `outf`: Open file to write to
    """
    for chr_name in coverage:
        outf.write('track type=wiggle_0 name=%s_PLUS description="%s PLUS" visibility=full color=0,0,255\n'%(
            title, description))
        outf.write("fixedStep chrom=%s start=1 step=1\n"%chr_name)
        for c in coverage[chr_name]['+']:
            outf.write("%d\n"%c)
        

        outf.write('track type=wiggle_0 name=%s_MINUS description="%s MINUS" visibility=full color=255,0,0\n'%(
            title, description))
        outf.write("fixedStep chrom=%s start=1 step=1\n"%chr_name)
        for c in coverage[chr_name]['-']:
            outf.write("%d\n"%c)


# Functions for mapping of chimeric fragments
def read_transcripts(trans_gff, feature='exon', identifier='gene_id'):
    """
    Read the transcripts, return a dictionary from transcript name to
    position
    Arguments:
    - `trans_gff`: A gff file of transcripts
    - `feature`: the name of the feature to index
    - `identifier`: The identifier to use from column 8
    """
    tus = {}
    for line in csv.reader(open(trans_gff), delimiter='\t'):
        if line[0].startswith('#') or line[0].startswith('@'):
            continue
        if line[2] != feature:
            continue
        ids_dict = {}
        for id_pair in line[8].strip().split(';'):
            try:
                k, v = id_pair.strip().split(' ')
            except ValueError:
                pass
            ids_dict[k] = v.replace('"','')
        fid = ids_dict[identifier]
        tus[fid] = (int(line[3])-1, int(line[4]), line[6])
    return tus


def get_unmapped_reads(
    samfile, outfile1, outfile2, length, maxG, rev=False, all_reads=False,
    dust_thr=0):
    """
    Get the list of unmapped paired reads and write the reads (mate 1 and 2) to
    the fastq files outfile1 and outfile2. The names of the reads is the same
    (assume equal in bam file)
    If rev is set assume first read is the reverse complement and reverse
    complement it, put it as read 2 and treat the second read as read 1.
    Can handle single-end as well.
    If all_reads is True, return the names of the reads that are mapped.
    Arguments:
    - `samfile`: Open Samfile object
    - `outfile1`: Open fastq file for reads 1
    - `outfile2`: Open fastq file for reads 2
    - `length`: Write the first X nt of the sequences
    - `maxG`: Maximal fraction of G's in any of the reads
    - `rev`: Reads are reverse complement (Livny's RNAtag-seq protocol for
             instance).
    - `all_reads`: Return all reads, including mapped ones
    - `dust_thr`: DUST filter threshold. If=0, not applied.
    """
    single_mapped = set()
    for read in samfile.fetch(until_eof=True):
        if (not read.is_paired) and (read.is_unmapped or all_reads):
            reverse_seq = False
            if read.is_reverse:
                # This can't happen unless all_reads is set to True
                reverse_seq = True
            cseq = read.seq
            cqual = read.qual
            # If the read is the reverse complement of the RNA XOR it's been
            # reversed on the bam file, reverse it
            if rev!=reverse_seq:
                cseq = str(Seq(cseq).reverse_complement())
                cqual = cqual[::-1]
            if all_reads and (not read.is_unmapped):
                single_mapped.add(read.qname)
            search_for = 'G'
            if rev:
                search_for = 'C'
            if cseq.count(search_for, 0, length) >= int(maxG*length) or\
                    cseq.count(search_for, -length) >= int(maxG*length):
                continue
            outfile1.write("@%s\n%s\n+\n%s\n"%(
                    read.qname, cseq[:length],
                    cqual[:length]))
            outfile2.write("@%s\n%s\n+\n%s\n"%(
                    read.qname, cseq[-length:],
                    cqual[-length:]))
            continue
        if (all_reads or read.is_unmapped or read.mate_is_unmapped or\
                (not read.is_proper_pair)) and read.is_paired:
            if all_reads and not (read.is_unmapped or read.mate_is_unmapped or\
                (not read.is_proper_pair)):
                single_mapped.add(read.qname)
            if read.is_read1==rev:
                ouf = outfile2
                outseq = Seq(read.seq)
                outqual = read.qual[-length:]
                # Reverse complement the read if it haven't been
                # done in the bam file. Otherwise, do nothing
                if not read.is_reverse:
                    outseq = outseq.reverse_complement()
                    outqual = read.qual[::-1][-length:]
                outseq = str(outseq[-length:])
                if (str(outseq).count('C')>=int(maxG*length)):
                    continue
            else: # First read in the fragment
                ouf = outfile1
                outseq = Seq(read.seq)
                outqual = read.qual[:length]
                if read.is_reverse:
                    outseq = outseq.reverse_complement()
                    outqual = read.qual[::-1][:length]
                outseq = str(outseq[:length])
                if outseq.count('G') >= int(maxG*length):
                    continue
            # test if read passes DUST filter
            if pass_dust_filter(outseq, dust_thr):
                ouf.write("@%s\n%s\n+\n%s\n"%(read.qname, outseq, outqual))
    return single_mapped

def pass_dust_filter(seq, thr):
    """
    Run dust filter and return True of False. The filter is applied as
    described in: ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/windowmasker_suppl.pdf
    The score is normalized to the length of the sequence (-3) and multiplied
    by 10 as the formula says 
    Arguments:
    - `seq`: The sequence
    - `thr`: The threshold
    """
    if len(seq) <= 3 or thr == 0:
        return True
    counts = defaultdict(int)
    for i in range(len(seq)-2):
        counts[seq[i:i+3]] += 1
    cscore = 0
    for c in counts.values():
        cscore += c * (c-1) * 0.5
#    sys.stderr.write("%g\n"%(cscore/(len(seq)-3)*10))
    return cscore/(len(seq)-3) * 10 <= thr

            
def run_dust_filter(fname, fout, prinseq_cmd, threshold):
    """
    Run the dust filter implemented in prinseq. Read a fastq file and write to
    the open file in fout.
    Arguments:
    - `fname`: Input file name
    - `fout`: Write output to this file
    - `prinseq_cmd`: prinseq.pl executable
    - `Threshold`: Above this value, reads will be ignored
    """
    call([prinseq_cmd, '-fastq', fname, '-lc_method', "dust", '-lc_threshold',
          str(threshold), '-out_good', 'stdout'], stdout=fout)
    


def get_XA_mapping(tags, max_mm=None):
    """
    Return a list of positions from the XA sam file data, these are additional
    mappings
    Arguments:
    - `tags`: The tags argument
    """
    tlist = {}
    for tpair in tags:
        if tpair[0] == 'XA':
            for xadat in tpair[1].split(';')[:-1]:
                alt_dat = xadat.split(',')
                if max_mm is None or int(alt_dat[-1])<=max_mm:
                    tlist[int(alt_dat[1][1:])] = alt_dat
                    
    return [tlist[i] for i in sorted(tlist)]


def get_NM_number(tags):
    """
    Return the number of mismatches found in the NM tag
    Arguments:
    - `tags`: The tags argument
    """
    for tpair in tags:
        if tpair[0] == 'NM':
            return tpair[1]
    return 10


def replace_with_XA(read, al, chrnames_bam):
    """
    Replace the position of the read with one of the alternatives
    Arguments:
    - `read`: A read object
    - `al`: Alternative position string
    - `chrnames_bam`: The names of the chrs in the bam file
    """
    apos = int(al[1][1:])
    nm_num = get_NM_number(read.tags)
    # Add the read one to the XA tag
    tags = read.tags
    for xt in tags:
        if xt[0] == 'XA':
            xaval = xt[1]
            tags.remove(xt)
            strs = '+'
            if read.is_reverse:
                strs = '-'
            tags.append(('XA', '%s,%s%d,%s,%d;'%(
                        chrnames_bam[read.tid],
                        strs, read.pos,
                        read.cigarstring, nm_num) +xaval))
            read.tags = tags
    read.pos = apos
    read.is_reverse = al[1][0]=='-'
    read.cigarstring = al[2]
    read.tid = min([i for i, x in enumerate(chrnames_bam) if x==al[0]])


def test_concordance(
    read1, read2, maxdist, chrnames_bam, trans_gff=None, remove_self=False):
    """
    Test if the two reads can be concordant and not ligated.
    If it finds a combination that is concordant, replace the concordant
    positions to be the major ones. You should check if they are under the
    allowed mismatches.
    Update:
    test all the XA, if one combination is concordant return True
    Arguments:
    - `read1`: read 1 object
    - `read2`: read 2 object
    - `maxdist`: Maximal distance to consider concordance if not on the same
                 transcript and have same orientation (possibly self ligation)
    - `chrnames_bam`: A list of chr names as they appear in the bam file.
                      can be generated from Samfile.getrname() function
    - `trans_gff`: A dict IGT->(from, to, strand)
    - `remove_self`: Remove pairs that have the same orientation but different
                     order i.e. r1--->r2---> instead of r2--->r1--->
    """
    def is_conc(str1, str2, pos1, pos2, chr1, chr2):
        """
        """
        if str1 != str2:
            return False
        if chr1 != chr2:
            return False
        if abs(pos2-pos1)<maxdist and ((pos1>pos2)==str1 or remove_self):
            return True
        if trans_gff:
            in_trans = False
            for tu_pos in trans_gff.values():
                if (tu_pos[2]=='-')==str1:
                    if tu_pos[0]<=pos1<=tu_pos[1] and tu_pos[0]<=pos2<=tu_pos[1]:
                        in_trans = True
                        break
            if in_trans and ((pos1>pos2)==str1 or remove_self):
                return True
        return False
                        
    
    # read1 is the reverse of the real read. If both directions are the same
    # and distance is short they can be concordant
    
    if read1.is_reverse == read2.is_reverse and read1.tid==read2.tid:
        if is_conc(
            read1.is_reverse, read2.is_reverse, read1.pos, read2.pos,
            read1.tid, read2.tid):
            return True
    r1_XA = get_XA_mapping(read1.tags)
    r2_XA = get_XA_mapping(read2.tags)
    if r1_XA:
        for altp in r1_XA:
            is_rev1 = (altp[1][0] == '-')
            pos1 = abs(int(altp[1]))
            if is_conc(
                is_rev1, read2.is_reverse, pos1, read2.pos, altp[0],
                chrnames_bam[read2.tid]):
                # Replace with alternative
                replace_with_XA(read1, altp, chrnames_bam)
                return True
            for altp2 in r2_XA:
                is_rev2 = (altp2[1][0] == '-')
                pos2 = abs(int(altp2[1]))
                if is_conc(is_rev1, is_rev2, pos1, pos2, altp[0], altp2[0]):
                    # Replace both with alternatives.
                    replace_with_XA(read1, altp, chrnames_bam)
                    replace_with_XA(read2, altp2, chrnames_bam)
                    return True
    if r2_XA:
        for altp2 in r2_XA:
            is_rev2 = (altp2[1][0] == '-')
            pos2 = abs(int(altp2[1]))
            if is_conc(read1.is_reverse, is_rev2, read1.pos, pos2, altp2[0], chrnames_bam[read1.tid]):
                # Replace with alternative
                replace_with_XA(read2, altp2, chrnames_bam)
                return True
        
    return False

def read_bam_file(bamfile, chrnames_bam, max_NM=0):
    """
    Given a Samfile object of a single mapped file, return the reads in the file
    If there are more than one mapping to a read return the first in the genome.
    Test if the read has less than (or equals to) number of allowed mismatches
    Arguments:
    - `bamfile`: A Samfile object
    - `chrnames_bam`: The names of the chromosomes
    - `max_NM`: Maximal number of mismatches
    """
    read_objects = {}
    for read in bamfile.fetch():
        if not read.is_unmapped:
            nm_num = get_NM_number(read.tags)
            if nm_num > max_NM:
                continue
            # If there are multiple hits, choose the one with the smallest coor
            alt_list = get_XA_mapping(read.tags, max_NM)
            min_pos = read.pos
            min_al = None
            for al in alt_list:
                apos = int(al[1][1:])
                # test this alternative only if its NM is as the original one
                if int(al[3])>nm_num:
                    continue
                if apos < min_pos:
                    min_pos = apos
                    min_al = al
            # If changed, add the read one to the XA tag
            if read.pos != min_pos:
                replace_with_XA(read, min_al, chrnames_bam)
            read_objects[read.qname] = read
    return read_objects


def write_reads_table(
    outfile, read1_reads, read2_reads, chrnames_bam, maxdist,
    remove_self, trans_gff=None, write_single=None, single_mapped=None,
    max_NM=1):
    """
    Read the lists of reads and print a list of chimeric fragments after
    removing concordant reads
    Arguments:
    - `outfile`: Print the reads positions to this open file
    - `read1_reads`: A dictionary of reads from side 1
    - `read2_reads`: A dictionary of reads from side 2, keys of 1 and 2 should
                     match
    - `chrnames_bam`: A list of chromosome names in the bam file
    - `maxdist`: Maximal distance between concordant reads
    - `remove_self`: Remove circular RNAs
    - `trans_gff`: A dictionary with transcripts positions, optional
    - `write_single`: Write reads that are not chimeric to this file
    - `single_mapped`: A set with read names of fragments that were originally
                       mapped as single
    - `max_NM`: Maximla number of mismatches. Used to screen printing of singles
    """
    

    for rname in read1_reads:
        if rname not in read2_reads:
            continue
        # If the two reads share at least one gene they are excluded
        write_to = outfile
        read_type = "chimera"
        if test_concordance(
            read1_reads[rname], read2_reads[rname], maxdist,
            chrnames_bam, trans_gff=trans_gff, remove_self=remove_self):
            if write_single:
                if get_NM_number(read1_reads[rname].tags) > max_NM or\
                        get_NM_number(read2_reads[rname].tags) > max_NM:
                    continue
                write_to = write_single
                read_type = "single"
            else:
                continue
        else:
            # If it was originally mapped as single and now as chimeric
            # ignore this read
            if single_mapped and rname in single_mapped:
                continue
        read1_chrn = chrnames_bam[read1_reads[rname].tid]
        read2_chrn = chrnames_bam[read2_reads[rname].tid]
        if not read2_reads[rname].is_reverse:
            end2_pos = read2_reads[rname].pos+read2_reads[rname].qlen-1
            end2_str = '+'
        else:
            end2_pos = read2_reads[rname].pos
            end2_str = '-'
        if read1_reads[rname].is_reverse:
            end1_pos = read1_reads[rname].pos+read1_reads[rname].qlen-1
            end1_str = '-'
        else:
            end1_pos = read1_reads[rname].pos
            end1_str = '+'
        write_to.write(
            "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n"%(
                read1_chrn, end1_pos+1, end1_str, read2_chrn, end2_pos+1,
                end2_str, rname, read_type))


def read_reads_table(reads_in, seglen, rRNAs=None, only_single=False):
    """
    Read a reads table and count the number of times each pair of segments
    appears in a chimeric fragment.
    Arguments:
    - `reads_in`: The table (tab delimited) of reads
    - `seglen`: The length of the segment
    - `rRNAs`: Remove chimeras that map to rRNA gene. This parameter holds the
               list of rRNA genes positions (chr, from, to ,strand)
    """
    region_interactions = defaultdict(lambda:defaultdict(list))
    region_ints_as1 = defaultdict(int)
    region_ints_as2 = defaultdict(int)
    total_interactions = 0
    for line in reads_in:
        end1_chrn, end1_pos1, end1_str, end2_chrn, end2_pos1, end2_str, _, rtype =\
            line.strip().split()
        end1_pos = int(end1_pos1)-1
        end2_pos = int(end2_pos1)-1
        if rRNAs:
            has_rRNA = False
            for rrgene in rRNAs:
                if end1_chrn == rrgene[0] and end1_str == rrgene[3] and\
                        rrgene[1] <= end1_pos < rrgene[2]:
                    has_rRNA = True
                    break
                if end2_chrn == rrgene[0] and end2_str == rrgene[3] and\
                        rrgene[1] <= end2_pos < rrgene[2]:
                    has_rRNA = True
                    break
            if has_rRNA:
                continue
        end1_seg = (end1_pos/seglen)*seglen
        end2_seg = (end2_pos/seglen)*seglen
        total_interactions += 1
        if (rtype != "single") != only_single:
            region_interactions[(end1_seg, end1_str, end1_chrn)]\
                [(end2_seg, end2_str, end2_chrn)].append((end1_pos, end2_pos))
        region_ints_as1[(end1_seg, end1_str, end1_chrn)] += 1
        region_ints_as2[(end2_seg, end2_str, end2_chrn)] += 1
    return (
        region_interactions, region_ints_as1, region_ints_as2,
        total_interactions)

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
        r1_from = int(line['Start of RNA1 first read'])-1
        r1_to = int(line['Start of RNA1 last read'])
        try:
            r1_chrn = chr_dict[line['RNA1 chromosome']]
        except KeyError:
            r1_chrn = line['RNA1 chromosome']
        r1_str = line['RNA1 strand']
        r2_from = int(line['Start of RNA2 last read'])-1
        r2_to = int(line['Start of RNA2 first read'])
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

def update_exp_in_vitro(
    region_interactions, region_ints_as1, region_ints_as2, sings_as_1,
    sings_as_2, cross_as_1, cross_as_2, ts_as_1, ts_as_2, vitro_ratio,
    factor_jump=10):
    """
    This method computes the expected number of in-vitro interactions using
    data from experiments with foreign RNA. Assuming the number of in-vitro
    interactions is proportional to the number of single interactions across
    experiments, one can divide the number of cross interactions of reg1 with
    the number of cross interactions of reg2, divide by the multiplication of
    singles of the two regions and multiply in the multiplication of the number
    of single reads in this library. The number of singles should be normalized
    to all single reads in the library. The number is then divided by
    cross_ratio which is the fraction of foreign RNA in the tube, if it's 0.5
    it means that the amount of self and foreign RNA was equal and if we have
    10 interactions with foreign, it might as well have another 10 with self.
    Arguments:
    - `region_interactions`: returned by read_reads_table
    - `region_ints_as1`: as above
    - `region_ints_as2`: as above
    - `sings_as_1`: returned from read_reads_table for the foreign singles
    - `sings_as_2`:as above
    - `cross_as_1`: reads for the cross ligation
    - `cross_as_2`: as above
    - `ts_as_1`: this library singles
    - `ts_as_2`: as above
    - `vitro_ratio`: A float with the estimate of in-vitro interactions out of
                     all interactions in the test library
                     RNA
    """
    total_counts = sum(region_ints_as1.values())
    all_counts = total_counts
    norm_factor_as_1 = {}
    norm_factor_as_2 = {}
    # The fragments as 1 can't also be as 2 since it's cross species
    sum_cross_as_1 = float(sum(cross_as_1.values()))
    sum_cross_as_2 = float(sum(cross_as_2.values()))
#    sum_all_cross = sum_cross_as_1 + sum_cross_as_2
    sum_singles_as_1 = float(sum(sings_as_1.values()))
    sum_singles_as_2 = float(sum(sings_as_2.values()))
    sum_ts_as_1 = float(sum(ts_as_1.values()))
    sum_ts_as_2 = float(sum(ts_as_2.values()))
    for region in cross_as_1:
        norm_factor_as_1[region] = ((cross_as_1[region]/sum_cross_as_1)/\
            ((sings_as_1[region] or 1)/sum_singles_as_1))*\
            (ts_as_1[region]/sum_ts_as_1)
    for region in cross_as_2:
        norm_factor_as_2[region] = ((cross_as_2[region]/sum_cross_as_2)/\
            ((sings_as_2[region] or 1)/sum_singles_as_2))*\
            (ts_as_2[region]/sum_ts_as_2)
    sum_of_vitro = 0
    vitro_counts = defaultdict(dict)
    mult_factor = vitro_ratio
    while sum_of_vitro < all_counts*min(vitro_ratio,1):
        sum_of_vitro = 0
        for reg1, r1data in region_interactions.items():
            if reg1 not in norm_factor_as_1:
                continue
            for reg2, its in r1data.items():
                if reg2 not in norm_factor_as_2:
                    continue
                counts = len(its)

                vitro_counts[reg1][reg2] = min(
                    int(norm_factor_as_1[reg1]*norm_factor_as_2[reg2]*all_counts*mult_factor),
                    counts)
                sum_of_vitro += vitro_counts[reg1][reg2]
        sys.stderr.write("%d\t%f\n"%(sum_of_vitro, mult_factor))
        mult_factor += factor_jump
    for reg1, r1data in vitro_counts.items():
        for reg2, ivt_counts in r1data.items():
            its_list = region_interactions[reg1][reg2]
            counts = len(its_list)
            counts -= ivt_counts
            region_interactions[reg1][reg2] = its_list[:counts]
            region_ints_as1[reg1] -= ivt_counts
            region_ints_as2[reg2] -= ivt_counts
            total_counts -= ivt_counts
    return total_counts



def minpv_regions(
    reg1, reg2, r_int, t_int_as1, t_int_as2, tot_int, f_int, seglen, maxsegs,
    min_odds):
    """
    return the regions that have minimal p-value. Exhaustive search
    Arguments:
    - `reg1`: Region1 seed
    - `reg2`: Region2 seed
    - `r_int`: regions interactions double dictionary
    - `t_int_as1`: Total interactions of every region as region 1
    - `t_int_as2`: As above for the second read
    - `tot_int`: Total interactions
    - `f_int`: interacting regions that were used
    - `seglen`: Length of segment (usually 100)
    - `maxsegs`: Maximal number of neighbouring segment to join
    - `min_odds`: Minimal odds-ratio to test
    """
    maxpv = 2
    maxpars = [0] * 9
    corr = 0
    for l1 in range(maxsegs):
        for l2 in range(maxsegs):
            for i1 in range(l1+1):
                for i2 in range(l2+1):
                    has_former = False
                    int_sum = 0
                    s1_sum = 0
                    s2_sum = 0
                    for r2 in range(l2+1):
                        real_r2 = reg2[0]+((r2-i2)*seglen)
                        s2_sum += t_int_as2[(real_r2, reg2[1], reg2[2])]
                    for r1 in range(l1+1):
                        real_r1 = reg1[0]+((r1-i1)*seglen)
                        s1_sum += t_int_as1[(real_r1, reg1[1], reg1[2])]
                        for r2 in range(l2+1):
                            real_r2 = reg2[0]+((r2-i2)*seglen)
                            if (real_r1, reg1[1], reg1[2], real_r2,reg2[1], reg2[2]) in f_int:
                                has_former = True
                                continue
                            int_sum += \
                                len(r_int[(real_r1, reg1[1], reg1[2])]\
                                        [(real_r2, reg2[1], reg2[2])])
                        
                    if has_former:
                        continue
                    # Set the values of the 2x2 contingency table a b c d:
                    #   a   b
                    #   c   d
                    a = int_sum
                    b = s1_sum - int_sum
                    c = s2_sum - int_sum
                    d = tot_int - s1_sum - s2_sum + int_sum
                    corr += 1
                    if b==0 or c==0:
                        odds = np.inf
                    else:
                        odds = (float(a)*d)/(float(b)*c)
                    if odds <= min_odds:
                        continue
                    odds, pv = fisher_exact(
                        [[a, b], [c, d]], alternative='greater')
                    if pv <= maxpv:
                        replace = False
                        if pv == maxpv:
                            # Both 0 probably, take the one with more
                            # interactions. If equal take the narrow one
                            if maxpars[4] < int_sum:
                                replace = True
                            else:
                                if (l1<=maxpars[0] and l2<maxpars[1]) or\
                                    (l2==maxpars[1] and l1<maxpars[0]):
                                    replace = True
                        else:
                            replace = True
                        if replace:
                            maxpv = pv 
                            maxpars = (l1, l2, i1, i2, int_sum, odds, b, c, d)
    return (maxpv*(corr or 1), maxpars[4], maxpars[5],
            reg1[0]-maxpars[2]*seglen,
            reg1[0]+(maxpars[0]-maxpars[2]+1)*seglen,
            reg2[0]-maxpars[3]*seglen,
            reg2[0]+(maxpars[1]-maxpars[3]+1)*seglen,
            maxpars[6], maxpars[7], maxpars[8])


# Functions for presenting the results
def read_targets(tarfile):
    """
    Return a tuple of targets
    Arguments:
    - `tarfile`: A tab-del file with EcoCyc names of sRNA and target
    """
    if not tarfile:
        return None
    tars = []
    try:
        for line in open(tarfile):
            tars.append(tuple(line.strip().split()[:2]))
    except IOError:
        return None
    return tars

def read_singles(singles_file):
    """
    Read the table of reads per gene and return a dictionary
    """
    if not singles_file:
        return None
    counts = {}
    try:
        for line in open(singles_file):
            spl = line.strip().split()
            try:
                counts[spl[0]] = sum([int(k) for k in spl[1:]])
            except ValueError:
                pass
    except IOError:
        return None
    return counts


def read_annotations(refseq_dir, an_ext = ('.ptt.gz', '.rnt.gz')):
    """
    Read the annotations from rnt and ptt files and return a dictionary
    Arguments:
    - `refseq_dir`: The Refseq dictionary
    - `an_ext`: Extensions of annotation files (iterable)
    """
    annotations = {}
    ec_files = []
    for ext in an_ext:
        ec_files.extend(glob.glob("%s/*%s"%(refseq_dir, ext)))
    for fin in ec_files:
        fo = gzip.open(fin)
        for row in csv.reader(fo, delimiter='\t'):
            try:
                annotations[row[4]] = row[8]
            except IndexError:
                pass
    return annotations


def get_genes_dict(ec_dir, pad=100):
    """
    Return a dictionary from genomic position to annotation based on EcoCyc
    Arguments:
    - `ec_dir`: EcoCyc directory
    - `pad`: Estimated UTR length
    """
    try:
        pos_maps, uid_gene = ecocyc_parser.get_mapping(ec_dir, pad)
    except IOError:
        return None
    pos_maps_lists = defaultdict(dict)
    for chrn, pos_data in pos_maps.items():
        for k, v in pos_data.items():
            try:
                cn1 = uid_gene[v[0]]['COMMON-NAME']
            except KeyError:
                cn1 = v[0]

            if len(v)>2 and (v[2]=='IGR' or v[2]=='IGT' or v[2]=='IGT_AS'):
                try:
                    cn2 = uid_gene[v[1]]['COMMON-NAME']
                except KeyError:
                    cn2 = ''
                outlist = list(v[:2])+[cn1, cn2, v[2]]
            else:
                outlist = list(v[:1])+[cn1]
                if len(v)>1:                                     
                    outlist.append(v[1])
            pos_maps_lists[chrn][k] = '.'.join(outlist)
    return pos_maps_lists


def list_of_genes(
    r1_region_from, r1_region_to, r1_str, r1_chrn, r1_chrnEC, r2_region_from,
    r2_region_to, r2_str, r2_chrn, r2_chrnEC, region_interactions, genes_dict,
    seglen, rlen):
    """
    Return a list of genes that intersect with the given regions.
    Return the minimal and maximal positions of reads in the regions as well.
    The genes are sorted according to the number of reads starting in each
    position
    Changed on Version 0.18:
    If a region overlaps antisense and another gene always prefer the other gene
    Arguments:
    - `r1_region_from`: the first position of r1 
    - `r1_region_to`: The end of r1
    - `r1_str`: r1 strand
    - `r1_chrn`: r1 chromosome
    - `r1_chrnEC`: r1 chromosome in EcoCyc mapping
    - `r2_region_from`: region 2 from
    - `r2_region_to`: region 2 to
    - `r2_str`: r2 strand
    - `r2_chrn`: r2 chromosome
    - `r2_chrnEC`: r2 chromosome in EcoCyc
    - `region_interactions`: A double dictionary region1->region2->[(p1, p2)]
    - `genes_dict`: dictionary of genes
    - `seglen`: length of segments
    - `rlen`: length of read
    """
    cdict_r1 = defaultdict(int)
    cdict_r2 = defaultdict(int)
    min1_pos = np.inf
    max1_pos = 0
    min2_pos = np.inf
    max2_pos = 0
    for r1_reg in range(r1_region_from, r1_region_to, seglen):
        for r2_reg in range(r2_region_from, r2_region_to, seglen):
            for (r1p, r2p) in region_interactions[(r1_reg, r1_str, r1_chrn)]\
                    [(r2_reg, r2_str, r2_chrn)]:
                min1_pos = min(min1_pos, r1p)
                max1_pos = max(max1_pos, r1p)
                min2_pos = min(min2_pos, r2p)
                max2_pos = max(max2_pos, r2p)
                if genes_dict:
                    if r1_chrnEC in genes_dict:
                        if r1_str == '-':
                            drange = range(max(r1p-rlen,0), r1p)
                        else:
                            drange = range(r1p, r1p+rlen)
                        for i in drange:
                            cdict_r1[genes_dict[r1_chrnEC][(i, r1_str)]] += 1
                    if r2_chrnEC in genes_dict:
                        if r2_str == '-':
                            drange = range(r2p, r2p+rlen)
                        else:
                            drange = range(max(r2p-rlen,0), r2p)
                        for i in drange:
                            cdict_r2[genes_dict[r2_chrnEC][(i, r2_str)]] += 1
    # reduce the AS counts to 1
    for k in cdict_r1:
        if k[-1] == 'AS':
            cdict_r1[k] = 1
    for k in cdict_r2:
        if k[-1] == 'AS':
            cdict_r2[k] = 1
    
    genes1_list = sorted(cdict_r1, key=cdict_r1.get, reverse=True)
    genes2_list = sorted(cdict_r2, key=cdict_r2.get, reverse=True)
    return genes1_list, genes2_list, min1_pos, min2_pos, max1_pos, max2_pos

def get_seqs(chrn, pfrom, pto, pstrand, fsa_seqs, shuffles=0):
    """
    Get the sequence for which coordinates are given. if shuffles>0 return
    a dictionary with the real sequence as the key 'real' and the shuffled
    with keys of numbers
    """
    pseq = fsa_seqs[chrn][pfrom:pto]
    if pstrand == '-':
        pseq = pseq.reverse_complement()
    if not pseq:
        pseq = 'AAA'
    if shuffles<0:
        shf_seqs = {'real': pseq}
        return shf_seqs
    if shuffles>0:
        shf_seqs = {'real': pseq}
        for i in range(shuffles):
            shf_seqs[i] = dinuc_shuffle.shuffle_difreq(str(pseq))
        return shf_seqs
    return pseq

def get_names(gname, uid_names, annotations):
    """
    return the names of the gene
    Arguments:
    - `gname`: a gene EcoCyc accession number
    """
    try:
        p0_name = uid_names[gname[0]]['COMMON-NAME']
    except KeyError:
        p0_name = gname[0]
        p0_desc = '-'
    else:
        try:
            p0_desc = annotations[p0_name]
        except KeyError:
            p0_desc = '-'
    if len(gname)>2 and (gname[2]=='IGR' or gname[2]=='IGT' or\
                             gname[2]=='IGT_AS'):
        try:
            cn2 = uid_names[gname[1]]['COMMON-NAME']
        except KeyError:
            cn2 = gname[1]
        p0_name += '.%s.%s'%(cn2, gname[2])
        try:
            p1_desc = annotations[cn2]
        except KeyError:
            p1_desc = '-'
        p0_desc += ' : %s'%p1_desc
        
    elif len(gname)>=2:
        p0_name += '.%s'%gname[1]
        
    return p0_name, p0_desc

def has_rep(chrn, rfrom, rto, rstrand, rep_pos):
    """
    Return the names of the REP elements in the region
    Arguments:
    - `chrn`: chromosome name
    - `rfrom`: region from
    - `rto`: region to
    - `rstrand`: region strand
    - `rep_pos`: REP positions dictionary
    """
    reps = []
    for repel, repp in rep_pos.items():
        if repp[0] == chrn and repp[3] == rstrand:
            if repp[1]<rto and repp[2]>=rfrom:
                reps.append(repel)
    return reps


def report_interactions(
    region_interactions, outfile, interacting_regions, seglen, ec_dir, genome, ec_chrs,
    refseq_dir, targets_file, rep_file,  single_counts, shuffles, RNAup_cmd,
    servers, rlen, est_utr_lens, pad_seqs, totRNA_count, ip_tot_norm=0,
    total_reads_IP=0, total_reads_total=0):
    """
    Report the interactions with additional data such as genes in region, if
    it's a known target, number of single fragments count, binding energy
    Arguments:
    - `region_interactions`: The double dictionary containing the reads
    - `outfile`: An open file to write the data to
    - `interacting_regions`: A list of interacting regions as tuples
    - `seglen`: Segment length
    - `ec_dir`: Directory of EcoCyc flatfiles
    - `genome`: Genome fasta file. Used if ec_dir is ommited
    - `ec_chrs`: A list of chromosome names to build a dictionary
    - `refseq_dir`: The RefSeq directory to get the gene descriptions from
    - `targets_file`: A file with sRNA-target in EcoCyc IDs
    - `rep_file`: A file containing the REP elements table
    - `single_counts`: A file with single counts, take the sum of each row
    - `shuffles`: Number of shuffles to compute empirical p-value
    - `RNAup_cmd`: command line of RNAup
    - `servers`: Number of CPUs or a list of servers to use.
    - `rlen`: Length of reads
    - `est_utr_lens`: Estimated lengths of UTRs when data is not available in
                      EcoCyc
    - `pad_seqs`: Pad the interacting regions when extracting sequences.
    - `totRNA_count`: A dictionary from the region (as in prev parameters) to
                      the number of reads from total RNA
    - `ip_tot_norm`: The maximal IP/total value. all values will be normalized
                     to this ratio
    - `total_reads_IP`: Number of reads in IP library
    - `total_reads_total`: Number of reads in total library
    """
    targets = read_targets(targets_file)
    singles = read_singles(single_counts)
    desc = read_annotations(refseq_dir)
#    genes_dict = get_genes_dict(ec_dir, est_utr_lens)
    try:
        pos_maps, _ = ecocyc_parser.get_mapping(ec_dir, est_utr_lens)
    except IOError:
        pos_maps = None
    try:
        rep_pos = ecocyc_parser.read_REP_table(rep_file)
    except IOError:
        rep_pos = None
    if genome and not ec_dir:
        fsa_seqs = {}
        from Bio import SeqIO
        for record in SeqIO.parse(genome, 'fasta'):
            fsa_seqs[record.id] = record.seq
    else: 
        try:
            fsa_seqs = ecocyc_parser.read_fsas(ec_dir)
        except IOError:
            fsa_seqs = None
    if shuffles != 0:
        rnup = RNAup_tar_pred.RNAupTarPred(cmd=RNAup_cmd, servers=servers)
    try:
        _, uid_names, _, sRNAs, _, _  = ecocyc_parser.read_genes_data(ec_dir)
    except IOError:
        uid_names, sRNAs = None, None
    if len(ec_chrs) >= 2 and ec_dir:
        chr_dict = dict(zip(ec_chrs.split(',')[0::2], ec_chrs.split(',')[1::2]))
    else:
        chr_dict = {}
    # All the output will be stored here and then sorted according to the name
    out_data = {}
    header_vec = [
        'RNA1 chromosome', 'Start of RNA1 first read', 'Start of RNA1 last read', 'RNA1 strand',
        'RNA2 chromosome', 'Start of RNA2 last read', 'Start of RNA2 first read', 'RNA2 strand',
        'interactions', 'other interactions of RNA1',
        'other interactions of RNA2', 'total other interactions', 'odds ratio',
        "Fisher's exact test p-value"]
    if ip_tot_norm > 0:
        header_vec.extend(
            ["total RNA reads1", "total RNA reads2", "lib norm IP RNA1",
             "lib norm IP RNA2", "lib norm total RNA1", "lib norm total RNA2",
             "IP/total ratio1", "IP/total ratio2",
             "Normalized Odds Ratio (NOR)",
             "RNA1 pred effect", "RNA2 pred effect", "Maximal RNA effect",
             "Total reads IP: %d"%total_reads_IP,
             "Total reads total: %d"%total_reads_total])
    if shuffles > 0 and fsa_seqs:
        header_vec.extend([
                'Free energy of hybridization',
                'empirical p-value of free energy'])
    if shuffles < 0 and fsa_seqs:
        header_vec.extend(['Free energy of hybridization'])
    if rep_pos:
        header_vec.extend(['RNA1 REP elements', 'RNA2 REP elements'])
    if targets and pos_maps:
        header_vec.append('Is known target')
    if singles and pos_maps:
        header_vec.extend([
                'RNA1 single fragments counts', 'RNA2 single fragments counts'])
    if pos_maps:
        header_vec = [
            'RNA1 EcoCyc ID', 'RNA2 EcoCyc ID', 'RNA1 name', 'RNA2 name',
            'RNA1 description', 'RNA2 description'] + header_vec
    # Used to sort the lines
    sorted_order = []
    sort_both_sRNAs = {}
    sort_one_sRNA = {}
    sort_3UTRs = {}
    sort_rest = {}

    for ittr in interacting_regions:
        (pv, ints, odds, r1_from, r1_to, r1_str, r1_chrnbam, r2_from, r2_to,
         r2_str, r2_chrnbam, mat_b, mat_c, mat_d) =\
            ittr
        rkey=(
            r1_from, r1_to, r1_str, r1_chrnbam, r2_from, r2_to, r2_str,
            r2_chrnbam)

        # get the genes in the region
        try:
            r1_chrn = chr_dict[r1_chrnbam]
        except KeyError:
            r1_chrn = r1_chrnbam
        try:
            r2_chrn = chr_dict[r2_chrnbam]
        except KeyError:
            r2_chrn = r2_chrnbam
        # Read genes and coordinates data
        genes1_list, genes2_list, min1_pos, min2_pos, max1_pos, max2_pos =\
            list_of_genes(
            r1_from, r1_to, r1_str, r1_chrnbam, r1_chrn, r2_from, r2_to, r2_str,
            r2_chrnbam, r2_chrn, region_interactions, pos_maps, seglen, rlen)
        out_data[rkey] = [
            r1_chrn, min1_pos+1, max1_pos+1, r1_str, r2_chrn, min2_pos+1,
            max2_pos+1, r2_str, ints, mat_b, mat_c, mat_d, odds, pv]

        if ip_tot_norm > 0:
            # Count the number of interactions of the regions with other regions
            tot_totals_as1 = sum([totRNA_count[(r1, r1_str, r1_chrnbam)] for\
                                      r1 in range(r1_from, r1_to, seglen)])
            tot_totals_as2 = sum([totRNA_count[(r2, r2_str, r2_chrnbam)] for\
                                      r2 in range(r2_from, r2_to, seglen)])
            rna1_eff = min(
                1,((mat_b + ints)/float(tot_totals_as1+1))/ip_tot_norm)*\
                ints/float(ints + mat_b)
            rna2_eff = min(
                1,((mat_c + ints)/float(tot_totals_as2+1))/ip_tot_norm) *\
                ints/float(ints + mat_c)
            pred_eff = min(
                1,((mat_b + ints)/float(tot_totals_as1+1))/ip_tot_norm)*\
                min(1,((mat_c + ints)/float(tot_totals_as2+1))/ip_tot_norm) *\
                odds
            out_data[rkey].extend(
                [tot_totals_as1, tot_totals_as2,
                 (ints+mat_b)/float(total_reads_IP),
                 (ints+mat_c)/float(total_reads_IP),
                 tot_totals_as1/float(total_reads_total),
                 tot_totals_as2/float(total_reads_total),
                 (mat_b + ints)/float(tot_totals_as1+1),
                 (mat_c + ints)/float(tot_totals_as2+1), pred_eff, rna1_eff,
                 rna2_eff, max(rna1_eff, rna2_eff)])
        if shuffles != 0 and fsa_seqs:
            p5_seqs = get_seqs(
                r1_chrn, min1_pos-pad_seqs, max1_pos+pad_seqs, r1_str, fsa_seqs,
                shuffles=shuffles)
            if r2_str == '+':
                p3_seq = get_seqs(
                    r2_chrn, min2_pos-pad_seqs, max2_pos, r2_str, fsa_seqs)
            else:
                p3_seq = get_seqs(
                    r2_chrn, min2_pos, max2_pos+pad_seqs, r2_str, fsa_seqs)
            rnrgs = rnup.scoreall(p3_seq, p5_seqs)
            if shuffles>0:
                pv = len([r for r in rnrgs.values() if r>=rnrgs['real']])/float(
                    len(rnrgs))
                out_data[rkey].extend([-rnrgs['real'], pv])
            else:
                out_data[rkey].extend([-rnrgs['real']])
        if rep_pos:
            out_data[rkey].append(','.join(has_rep(
                        r1_chrn, min1_pos, max1_pos, r1_str, rep_pos)))
            out_data[rkey].append(','.join(has_rep(
                        r2_chrn, min2_pos, max2_pos, r2_str, rep_pos)))
            
        if not pos_maps:
            sorted_order.append(rkey)
            continue
        if not genes1_list:
            genes1_list = ['-']
        if not genes2_list:
            genes2_list = ['-']
        gname1 = genes1_list[0]
        gname2 = genes2_list[0]
        g1common, g1desc = get_names(gname1, uid_names, desc)
        g2common, g2desc = get_names(gname2, uid_names, desc)
        out_data[rkey] = [
            '.'.join([str(j) for j in gname1]),
            '.'.join([str(j) for j in gname2]), g1common,
            g2common, g1desc, g2desc] + out_data[rkey]
        if targets:
            is_target = False
            for g1 in gname1:
                for g2 in gname2:
                    is_target |= (g1, g2) in targets or (g2, g1) in targets
            out_data[rkey].append(is_target)
        if singles:
            for gg in (gname1, gname2):
                gnums = []
                for g1 in gg:
                    try:
                        gnums.append(str(singles[g1]))
                    except KeyError:
                        pass
                out_data[rkey].append(':'.join(gnums))
            
        if gname1[0] in sRNAs and gname2[0] in sRNAs:
            sort_both_sRNAs[rkey] = ''.join(sorted([g1common, g2common]))
        elif gname1[0] in sRNAs:
            sort_one_sRNA[rkey] = g1common+g2common
        elif gname2[0] in sRNAs:
            sort_one_sRNA[rkey] = g2common+g1common
        elif '3UTR' in gname1 or 'EST3UTR' in gname1:
            sort_3UTRs[rkey] = g1common+g2common
        elif '3UTR' in gname2 or 'EST3UTR' in gname2:
            sort_3UTRs[rkey] = g2common+g1common
        else:
            sort_rest[rkey] = ''.join(sorted([g1common, g2common]))

    # Sort the list of interactions (if pos_maps is set)
    sorted_order.extend(sorted(sort_both_sRNAs, key=sort_both_sRNAs.get))
    sorted_order.extend(sorted(sort_one_sRNA, key=sort_one_sRNA.get))
    sorted_order.extend(sorted(sort_3UTRs, key=sort_3UTRs.get))
    sorted_order.extend(sorted(sort_rest, key=sort_rest.get))
    # Print all the interactions
    outer=csv.writer(outfile, delimiter='\t')
    outer.writerow(header_vec)
    for k in sorted_order:
        outer.writerow(out_data[k])
