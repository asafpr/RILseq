"""
Parse the genes.dat, promoters.dat, terminators.dat files to generate data on
coverage of the genome and genes names. In addition has a function for reading
the genomic sequences.
"""

from collections import defaultdict
import csv

def get_mapping(ec_dir='/home/users/assafp/Database/EcoCyc/current/data',
                utr_len=100, inc_rep_pos=False):
    """
    Get the mapping directly
    Arguments:
    - `ec_dir`: EcoCyc dir
    - `utr_len`: For estimating UTRs
    - `inc_rep_pos`: Include REP positions in map?
    """
    tu_promoters = read_promoters_data(ec_dir)
    tu_terminators = read_terminators_data(ec_dir)
    if inc_rep_pos:
        rep_pos = read_REP_table(ec_dir)
    else:
        rep_pos = None
    uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list, rRNAs =\
        read_genes_data(ec_dir)
    tu_genes = defaultdict(list)
    for gene, tus in uid_tudata.items():
        for tu in tus:
            tu_genes[tu].append(gene)
    fsa_seqs = read_fsas(ec_dir)
    fsa_lens = {}
    for fn, fs in fsa_seqs.items():
        fsa_lens[fn] = len(fs)
    maps = position_to_gene(
        uid_pos, uid_tudata, tu_promoters, tu_terminators, tu_genes,
        sRNAs_list+other_RNAs_list, fsa_lens, rep_pos=rep_pos, utr_len=utr_len)
    return maps, uid_names

def read_REP_table(rep_file):
    """
    Read the REP elements positions from a txt file. This file was generated
    using SmartTables
    (e.g. this: http://ecocyc.org/group?id=biocyc14-8223-3640227683)
    Arguments:
    - `rep_file`: The table files, coordinates are in columns 3,4,6,7
    """
    rep_pos = {}
    for line in csv.reader(open("%s"%rep_file), delimiter='\t'):
        chrn = line[6].rsplit('-',1)[0]
        try:
            rep_pos[line[0]] = [chrn, int(line[2])-1, int(line[3]), line[5]]
        except ValueError:
            pass
    return rep_pos
    
    


def read_fsas(ec_dir='/home/users/assafp/Database/EcoCyc/current/data'):
    """
    Read the fasta files of the chromosomes, return a dictionary from the
    chrname to the sequence
    Arguments:
    - `ec_dir`: main directory of EcoCyc
    """
    import glob
    from Bio import SeqIO
    fsa_seqs = {}
    for fsa_name in glob.glob("%s/*-*.fsa"%(ec_dir)):
        fsa_seqs[fsa_name.rsplit('/',1)[1].rsplit('.',1)[0].split('-',1)[1]] =\
                     SeqIO.read(open(fsa_name), 'fasta').seq
    return fsa_seqs

def generate_transcripts_file(
    outfile, utr_len=100,
    ec_dir='/home/users/assafp/Database/EcoCyc/current/data', chr_dict=None):
    """
    Write a gff file with transcripts boundaries
    Arguments:
    - `outfile`: An open out file
    - `utr_len`: Default UTR len
    - `ec_dir`: EcoCyc data dir
    - `chr_dict`: map chromosome name in EcoCyc to another name (e.g. chr)
    """
    tu_promoters = read_promoters_data(ec_dir)
    tu_terminators = read_terminators_data(ec_dir)
    uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list, rRNAs =\
        read_genes_data(ec_dir)
    tu_genes = defaultdict(list)
    for gene, tus in uid_tudata.items():
        if not tus:
            tu_genes[gene] = [gene]
        for tu in tus:
            tu_genes[tu].append(gene)
    tu_boundaries = {}
    for tu in tu_genes:
        tu_str = uid_pos[tu_genes[tu][0]][3]
        tu_chrn = uid_pos[tu_genes[tu][0]][0]
        if chr_dict:
            try:
                tu_chrn = chr_dict[tu_chrn]
            except KeyError:
                pass
        tu_boundaries[tu] = [0, 0, tu_str, tu_chrn]
        if tu in tu_promoters:
            if tu_str == '+':
                tu_boundaries[tu][0] = tu_promoters[tu]
            else:
                tu_boundaries[tu][1] = tu_promoters[tu]
        else:
            # Get first gene
            if tu_str == '+':
                first_pos = min([uid_pos[gene][1] for gene in tu_genes[tu]])
                tu_boundaries[tu][0] = first_pos-utr_len
            else:
                last_pos = max([uid_pos[gene][2] for gene in tu_genes[tu]])
                tu_boundaries[tu][1] = last_pos+utr_len
        if tu in tu_terminators:
            if tu_str == '+':
                tu_boundaries[tu][1] = max(tu_terminators[tu])+1
            else:
                tu_boundaries[tu][0] = min(tu_terminators[tu])
        else:
            # Get first gene
            if tu_str == '-':
                first_pos = min([uid_pos[gene][1] for gene in tu_genes[tu]])
                tu_boundaries[tu][0] = first_pos-utr_len
            else:
                last_pos = max([uid_pos[gene][2] for gene in tu_genes[tu]])
                tu_boundaries[tu][1] = last_pos+utr_len
    # Now write all the tus to a gff file
    for tu, tub in tu_boundaries.items():
        outfile.write(
            '%s\tEcoCyc\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(tub[3], tub[0]+1, tub[1], tub[2], tu, tu))
        

def generate_gff_file(
    outfile, ec_dir='/home/users/assafp/Database/EcoCyc/current/data',
    chr_dict=None):
    """
    Generate a gff file using the Unique ID
    Arguments:
    - `outfile`: An open file
    - `ec_dir`: EcoCyc data dir
    - `chr_dict`: map chromosome name in EcoCyc to another name (e.g. chr)
    """
    uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list, rRNAs =\
        read_genes_data(ec_dir)
    sort_dict = {}
    for gi, gpos in uid_pos.items():
        sort_dict[gi] = gpos[1]
    for gname in sorted(sort_dict, key=sort_dict.get):
        gpos = uid_pos[gname]
        if chr_dict:
            try:
                chrn = chr_dict[gpos[0]]
            except KeyError:
                chrn = gpos[0]
        outfile.write(
            '%s\tEcoCyc\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n'%(chrn, gpos[1]+1, gpos[2], gpos[3], gname, gname))
         

def read_promoters_data(
    ec_dir='/home/users/assafp/Database/EcoCyc/current/data',
    pfile = 'promoters.dat'):
    """
    Read the promoters dat file and return a dictionary tu-> promoter +1
    position
    Arguments:
    - `ec_dir`: The EcoCyc data directory
    - `pfile`: The promoters data file name
    """
    p1pos = None
    tuname = None
    tu_promoters = {}
    with open("%s/%s"%(ec_dir, pfile)) as pin:
        for line in pin:
            if line.startswith('#'):
                continue
            if line.startswith('ABSOLUTE-PLUS-1-POS'):
                p1pos = int(line.strip().split()[-1]) - 1
            if line.startswith('COMPONENT-OF'):
                tutemp = line.strip().split()[-1]
                if tutemp.startswith('TU'):
                    tuname = tutemp.replace('TU', 'IGT')
            if line.startswith('//'):
                if tuname and p1pos:
                    tu_promoters[tuname] = p1pos
                p1pos = None
                tuname = None
    return tu_promoters


def read_terminators_data(
    ec_dir='/home/users/assafp/Database/EcoCyc/current/data',
    tfile = 'terminators.dat'):
    """
    Read the promoters dat file and return a dictionary
    tu-> (term left, term right)
    position
    Arguments:
    - `ec_dir`: The EcoCyc data directory
    - `tfile`: The terminators data file name
    """
    left_pos = None
    right_pos = None
    tuname = None
    tu_promoters = {}
    with open("%s/%s"%(ec_dir, tfile)) as tin:
        for line in tin:
            if line.startswith('#'):
                continue
            if line.startswith('LEFT-END-POSITION'):
                left_pos = int(line.strip().split()[-1]) - 1
            if line.startswith('RIGHT-END-POSITION'):
                right_pos = int(line.strip().split()[-1]) - 1
            if line.startswith('COMPONENT-OF'):
                tutemp = line.strip().split()[-1]
                if tutemp.startswith('TU'):
                    tuname = tutemp.replace('TU', 'IGT')
            if line.startswith('//'):
                if tuname and left_pos and right_pos:
                    tu_promoters[tuname] = (left_pos, right_pos)
                left_pos = None
                right_pos = None
                tuname = None
    return tu_promoters



def read_genes_data(
    ec_dir='/home/users/assafp/Database/EcoCyc/current/data',
    gfile = 'genes.dat', sRNAs_types=('BC-2.2','BC-2.2.7'),
    other_RNAs_types=('BC-2.2.5', 'BC-2.2.6', 'BC-3.1.3.6'),
    exclude_set=('EG10438',),
    rRNA_prod='RRNA'):
    """
    Read the genes.dat from the ec directory 
    Arguments:
    - `ec_dir`: The EcoCyc data directory
    - `gfile`: The genes data file name
    - `sRNAs_type`: return a list of genes with this TYPES identifier
    - `other_RNAs_type`: treat these genes as RNA genes (for UTRs)
    - `exclude_set`: Exclude these genes from the sRNAs list (hfq)
    Return:
    - `uid_pos`: A dictionary UID->genomic position (from, to, strand) 0-based
                 The right position is exclusive as python coords
    - `uid_names`: A dict UID->{category->name}
    - `uid_tudata`: A dict UID->[IGTs]
    - `sRNAs list`: A list of sRNAs,
    - `other_RNAs_list`: a list of tRNAs, rRNAs etc
    - `rRNA_prod`: Name of rRNA product to return a list of rRNAs
    """
    uid_pos = {}
    uid_names = {}
    sRNAs_list = []
    other_RNAs_list = []
    rRNAs = []
    uid_tudata = {}
    with open("%s/%s"%(ec_dir, gfile)) as ecin:
        ingene = None
        has_types = False
        has_other_types = False
        accs = {}
        coords = [0,0,'']
        tu_data = []
        chrom_name = None
        for line in ecin:
            if line.startswith('#'):
                continue
            if line.startswith('UNIQUE-ID'):
                ingene = line.strip().split()[-1]
            if line.startswith('TYPES'):
                has_types |= line.strip().split()[-1] in sRNAs_types
                has_other_types |= line.strip().split()[-1] in other_RNAs_types 
            if line.startswith('COMMON-NAME') or line.startswith('ACCESSION'):
                spl = line.strip().split()
                accs[spl[0]] = spl[2]
            if line.startswith('LEFT-END-POSITION'):
                coords[0] = int(line.strip().split()[-1]) - 1
            if line.startswith('RIGHT-END-POSITION'):
                coords[1] = int(line.strip().split()[-1])
            if line.startswith('TRANSCRIPTION-DIRECTION'):
                coords[2] = line.strip().split()[-1]
            if line.startswith('COMPONENT-OF'):
                tuname = line.strip().split()[-1]
                if tuname.startswith('TU'):
                    tu_data.append(tuname.replace('TU', 'IGT'))
                else:
                    chrom_name = tuname.rsplit('-',1)[0]
            if line.startswith('PRODUCT'):
                if line.strip().split()[-1].split('-')[-1]==rRNA_prod:
                    rRNAs.append(ingene)
            if line.startswith('//'):
                if coords[2] != '' and coords[1]>coords[0]:
                    uid_pos[ingene] = [chrom_name] + coords
                    uid_names[ingene] = accs
                    uid_tudata[ingene] = set(tu_data)
                    if has_types and (ingene not in exclude_set): #exclude hfq
                        sRNAs_list.append(ingene)
                    if has_other_types and (ingene not in exclude_set):
                        other_RNAs_list.append(ingene)
                ingene = None
                has_types = False
                has_other_types = False
                accs = {}
                coords = [0,0,'']
                tu_data = []
                chrom_name = None
                
    return  uid_pos, uid_names, uid_tudata, sRNAs_list, other_RNAs_list, rRNAs

def anti_strand(strand):
    """
    Return the opposite strand
    Arguments:
    - `strand`: + or -
    """
    if strand == '+':
        return '-'
    else:
        return '+'


def position_to_gene(
    uid_pos, uid_tu, tu_promoters, tu_terminators, tu_genes, sRNAs_list,
    fsa_lens, rep_pos=None, utr_len=100):
    """
    For each postion (and strand) determines the gene in the position.
    If the real transcription is known use it. If between the TSS and the
    gene start there is another gene, both are in the same IGT,
    this region is intergenic and will be marked as gene1_gene2_IGT.
    If the terminator is known use it to determine 3' end of the IGT.
    If the gene is a sRNA use the end of the gene as 3' end
    If there are multiple IGTs associated with the gene take the longest one.
    If the promoter/terminator is not known (not associated in the DB)
    use the default length to set 3'/5' UTRs.
    If the region is intergenic annotate is as gene1_gene2_IGR
    Arguments:
    - `uid_pos`: Dict Unique ID -> position
    - `uid_tu`: Unique ID -> Transcription units
    - `tu_promoters`: IGT ID -> promoter position (+1)
    - `tu_terminators`: IGT ID -> terminator positions (left, right)
    - `tu_genes`: Dict IGT ID -> [genes]
    - `sRNAs_list`: A list of sRNAs, ignore terminators and promoters
    - `fsa_lens`: Lengths of chromosomes
    - `rep_pos`: Positions of REP elements in the genome
    - `utr_len`: Default UTR length

    Returns:
    - `pos_map`: A dictionary (position, strand) -> annotation
    """
    def first_in_TU(gname):
        """
        If it's first in TU return the promoter name (IGT). If not return None
        Arguments:
        - `gname`: gene name
        """
        strand = uid_pos[gname][3]
        chrn = uid_pos[gname][0]
        longest = fsa_lens[chrn]
        if strand == '-':
            longest = 0
        longest_name = None
        # Get the farthest promoter
        for tu in uid_tu[gname]:
            if tu in tu_promoters:
                if strand == '+':
                    if tu_promoters[tu] < longest:
                        longest = tu_promoters[tu]
                        longest_name = tu
                else:
                    if tu_promoters[tu] > longest:
                        longest = tu_promoters[tu]
                        longest_name = tu
        if not longest_name:
            return None
        # Test if there are genes between the gene and the promoer
        for alt_gene in tu_genes[longest_name]:
            if alt_gene == gname:
                continue
            gpos = uid_pos[alt_gene]
            if gpos[3] == strand:
                if strand == '+' and longest < gpos[1] < uid_pos[gname][1]:
                    return None
                elif uid_pos[gname][2] < gpos[2] < longest:
                    return None
        return longest_name
    
    def last_in_TU(gname):
        """
        If it's last in TU return the terminator name (IGT). If not return None
        Arguments:
        - `gname`: gene name
        """
        strand = uid_pos[gname][3]
        chrn = uid_pos[gname][0]
        longest = 0
        if strand == '-':
            longest = fsa_lens[chrn]
        longest_name = None
        # Get the farthest promoter
        for tu in uid_tu[gname]:
            if tu in tu_terminators:
                if strand == '+':
                    if max(tu_terminators[tu]) > longest:
                        longest = max(tu_terminators[tu])
                        longest_name = tu
                else:
                    if min(tu_terminators[tu]) < longest:
                        longest = min(tu_terminators[tu])
                        longest_name = tu
        if not longest_name:
            return None
        # Test if there are genes between the gene and the promoer
        for alt_gene in tu_genes[longest_name]:
            if alt_gene == gname:
                continue
            gpos = uid_pos[alt_gene]
            if gpos[3] == strand:
                if strand == '+' and longest > gpos[1] > uid_pos[gname][1]:
                    return None
                elif uid_pos[gname][2] > gpos[2] > longest:
                    return None
        return longest_name
    
    # First set all the annotations of genes (CDS)
    pos_map = defaultdict(dict)
    for gene, positions  in uid_pos.items():
        for i in range(positions[1], positions[2]):
            pos_map[positions[0]][(i, positions[3])] = (gene,)
            if (i, anti_strand(positions[3])) not in pos_map[positions[0]]:
                pos_map[positions[0]][(i, anti_strand(positions[3]))] =\
                    (gene, 'AS')
    # That was the easy part. Now move to the 5' UTRs
    for gname, gpos in uid_pos.items():
        if gname in sRNAs_list:
            continue
        # If no TU or no promoter or terminator is assigned to any of the TU
        # use the default length
        est_5utr = gname not in uid_tu
        est_3utr = gname not in uid_tu
        if not est_5utr:
            is_first = True
            has_utr = False
            for tu in uid_tu[gname]:
                has_utr |= tu in tu_promoters
                for alt_gene in tu_genes[tu]:
                    if gpos[3] == '+' and uid_pos[alt_gene][1] < gpos[1]:
                        is_first = False
                    elif gpos[3] == '-' and uid_pos[alt_gene][2] > gpos[2]:
                        is_first = False
            est_5utr = is_first and not has_utr
        if not est_3utr:
            is_last = True
            has_utr = False
            for tu in uid_tu[gname]:
                has_utr |= tu in tu_terminators
                for alt_gene in tu_genes[tu]:
                    if gpos[3] == '+' and uid_pos[alt_gene][2] > gpos[2]:
                        is_last = False
                    elif gpos[3] == '-' and uid_pos[alt_gene][1] < gpos[1]:
                        is_last = False
            est_3utr = is_last and not has_utr
        # Estimate 5' UTR
        if est_5utr:
            if gpos[3] == '+':
                for i in range(gpos[1]-utr_len, gpos[1]):
                    if (i, gpos[3]) not in pos_map[gpos[0]]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', gpos[1]-i)
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST5UTR_AS', gpos[1]-i)
                    # if it is, test if it's estimated and the length
                    # of this UTR is smaller than the other, already written
                    # UTR, if it's indeed shorter, replace it
                    elif len(pos_map[gpos[0]][(i, gpos[3])])>2 and \
                            pos_map[gpos[0]][(i, gpos[3])][1].startswith('EST') and \
                            pos_map[gpos[0]][(i, gpos[3])][2] > gpos[1]-i:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', gpos[1]-i)
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST5UTR_AS', gpos[1]-i)
                    elif pos_map[gpos[0]][(i, gpos[3])][-1] == 'AS':
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', gpos[1]-i)
            else: # Minus strand, do the same
                for i in range(gpos[2], gpos[2]+utr_len):
                    if (i, gpos[3]) not in pos_map[gpos[0]]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', i-gpos[2])
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST5UTR_AS', i-gpos[2])
                    elif len(pos_map[gpos[0]][(i, gpos[3])])>2 and\
                            pos_map[gpos[0]][(i, gpos[3])][1].startswith('EST') and\
                            pos_map[gpos[0]][(i, gpos[3])][2] > i-gpos[2]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', i-gpos[2])
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST5UTR_AS', i-gpos[2])
                    elif pos_map[gpos[0]][(i, gpos[3])][-1] == 'AS':
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST5UTR', i-gpos[2])

        if est_3utr:
            curr_chr_len = fsa_lens[gpos[0]]
            if gpos[3] == '+':
                for i in range(gpos[2], gpos[2]+utr_len):

                    i = i % curr_chr_len
                    if (i, gpos[3]) not in pos_map[gpos[0]]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', i-gpos[2])
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST3UTR_AS', i-gpos[2])
                    elif len(pos_map[gpos[0]][(i, gpos[3])])>2 and \
                            pos_map[gpos[0]][(i, gpos[3])][1].startswith('EST') and \
                            pos_map[gpos[0]][(i, gpos[3])][2] > i-gpos[2]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', i-gpos[2])
#                        if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                            pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                                gname, 'EST3UTR_AS', i-gpos[2])
                    elif pos_map[gpos[0]][(i, gpos[3])][-1] == 'AS':
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', i-gpos[2])

            else: # Minus strand
                for i in range(gpos[1]-utr_len, gpos[1]):
                    i = i % curr_chr_len
                    
                    if (i, gpos[3]) not in pos_map[gpos[0]]:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', gpos[1]-i)
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                            gname, 'EST3UTR_AS', gpos[1]-i)
                    elif len(pos_map[gpos[0]][(i, gpos[3])])>2 and \
                            pos_map[gpos[0]][(i, gpos[3])][1].startswith('EST') and \
                            pos_map[gpos[0]][(i, gpos[3])][2] > gpos[1]-i:
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', gpos[1]-i)
#                        if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                            pos_map[gpos[0]][(i, anti_strand(gpos[3]))] = (
#                                gname, 'EST3UTR_AS', gpos[1]-i)
                    elif pos_map[gpos[0]][(i, gpos[3])][-1] == 'AS':
                        pos_map[gpos[0]][(i, gpos[3])] =\
                            (gname, 'EST3UTR', gpos[1]-i)

        if est_5utr and est_3utr:
            continue
        # If there is a TU, map to the real promoter and terminator
        # Set the 5' UTR of gene
        if not est_5utr:
            first_tu_name = first_in_TU(gname)
        else:
            first_tu_name = None
        if first_tu_name:
            if gpos[3] == '+':
                for i in range(tu_promoters[first_tu_name], gpos[1]):
                    if (i, gpos[3]) not in pos_map[gpos[0]] or\
                            pos_map[gpos[0]][(i, gpos[3])][-1]=='AS':
                        pos_map[gpos[0]][(i, gpos[3])] = (gname, '5UTR')
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] =\
#                            (gname, '5UTR_AS')
            else:
                for i in range(gpos[2], tu_promoters[first_tu_name]+1):
                    if (i, gpos[3]) not in pos_map[gpos[0]] or\
                            pos_map[gpos[0]][(i, gpos[3])][-1]=='AS':
                        pos_map[gpos[0]][(i, gpos[3])] = (gname, '5UTR')
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] =\
#                            (gname, '5UTR_AS')
        # Set the 3' end of the gene
        if not est_3utr:
            last_tu_name = last_in_TU(gname)
        else:
            last_tu_name = None
        if last_tu_name:
            if gpos[3] == '+':
                for i in range(gpos[2], max(tu_terminators[last_tu_name])+1):
                    if (i, gpos[3]) not in pos_map[gpos[0]] or\
                            pos_map[gpos[0]][(i, gpos[3])][-1]=='AS':
                        pos_map[gpos[0]][(i, gpos[3])] = (gname, '3UTR')
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] =\
#                            (gname, '3UTR_AS')
            else:
                for i in range(min(tu_terminators[last_tu_name]), gpos[1]):
                    if (i, gpos[3]) not in pos_map[gpos[0]] or\
                            pos_map[gpos[0]][(i, gpos[3])][-1]=='AS':
                        pos_map[gpos[0]][(i, gpos[3])] = (gname, '3UTR')
#                    if (i, anti_strand(gpos[3])) not in pos_map[gpos[0]]:
#                        pos_map[gpos[0]][(i, anti_strand(gpos[3]))] =\
#                            (gname, '3UTR_AS')
        
    # Close the gaps with IGR regions
    for chrn, chrlen in fsa_lens.items():
        try:
            max_plus_pos = max([j for j in range(chrlen) if\
                                    (j, '+') in pos_map[chrn] and\
                                    len(pos_map[chrn][(j, '+')]) == 1])
        except ValueError:
            max_plus_pos = None
        try:
            max_min_pos = max([j for j in range(chrlen) if\
                                   (j, '-') in pos_map[chrn] and\
                                   len(pos_map[chrn][(j, '-')]) == 1])
        except ValueError:
            max_min_pos = None
        if max_plus_pos and max_plus_pos > max_min_pos:
            last_gene = pos_map[chrn][(max_plus_pos, '+')]
        elif max_min_pos:
            last_gene = pos_map[chrn][(max_min_pos, '-')]
        else:
            last_gene = None
        i = 0
#    import sys
#    sys.stderr.write('%s\n'%last_gene[0])
        while last_gene and i < chrlen:
            if (i, '+') in pos_map[chrn] and len(pos_map[chrn][(i, '+')]) == 1:
                last_gene = pos_map[chrn][(i, '+')]
            elif (i, '-') in pos_map[chrn] and len(pos_map[chrn][(i, '-')])== 1:
                last_gene = pos_map[chrn][(i, '-')]
            if (i, '+') in pos_map[chrn] and (i, '-') in pos_map[chrn]:
                i += 1
                continue
            next_pos = i+1
            next_strand = None
            while next_pos < chrlen:
                if (next_pos, '+') in pos_map[chrn] and\
                        len(pos_map[chrn][(next_pos, '+')])==1:
                    next_strand = '+'
                    break
                if (next_pos, '-') in pos_map[chrn] and\
                        len(pos_map[chrn][(next_pos, '-')])==1:
                    next_strand = '-'
                    break
                next_pos += 1
            if next_strand == '+':
                next_gene = pos_map[chrn][(next_pos, '+')]
            elif next_strand == '-':
                next_gene = pos_map[chrn][(next_pos, '-')]
            else:
                next_gene = ('None', )
#        sys.stderr.write('%d\t%s\n'%(i, next_gene[0]))
            igr_plus_tag = 'IGR'
            igr_minus_tag = 'IGR'
        # Test if the genes are in the same IGT
            for tu in uid_tu[last_gene[0]]:
                if next_gene[0] in tu_genes[tu]:
                    if uid_pos[last_gene[0]][3] == '+':
                        igr_plus_tag = 'IGT'
                        igr_minus_tag = 'IGT_AS'
                    else:
                        igr_plus_tag = 'IGT_AS'
                        igr_minus_tag = 'IGT'

            while (((i, '+') not in pos_map[chrn]) or\
                       ((i, '-') not in pos_map[chrn])) and i < chrlen:
                if (i, '+') not in pos_map[chrn]:
                    pos_map[chrn][(i, '+')] =\
                        (last_gene[0], next_gene[0], igr_plus_tag)
                if (i, '-') not in pos_map[chrn]:
                    pos_map[chrn][(i, '-')] =\
                        (last_gene[0], next_gene[0], igr_minus_tag)
                i += 1
        # Decorate with REP elements positions
        if rep_pos:
            for i in range(chrlen):
                for rep_name, rep_d in rep_pos.items():
                    if rep_d[0] == chrn and rep_d[1]<=i<rep_d[2]:
                        pos_map[chrn][(i, rep_d[3])] = tuple(
                            list(pos_map[chrn][(i, rep_d[3])]) +\
                                ['REP_%s'%rep_name])
                        pos_map[chrn][(i, anti_strand(rep_d[3]))] = tuple(
                            list(pos_map[chrn][(i, anti_strand(rep_d[3]))]) +\
                                ['REP_%s_AS'%rep_name])

    return pos_map

