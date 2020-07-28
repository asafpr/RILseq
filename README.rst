
.. image:: https://badge.fury.io/py/RILseq.svg
    :target: https://badge.fury.io/py/RILseq

.. image:: https://travis-ci.org/asafpr/RILseq.svg
    :target: https://travis-ci.org/asafpr/RILseq

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :target: http://bioconda.github.io/recipes/rilseq/README.html

================
RILseq
================
Intention
---------
This package can be used to analyzed RILseq experiments. It is written for a prokaryotic genome, without splice junction mapping and with some additional features. RILseq is described in Melamed et al, Molecular Cell 63 (2016), pp. 884-897 (http://www.cell.com/molecular-cell/fulltext/S1097-2765(16)30413-0) and Melamed et al, Nature Protocols volume 13, pages 133 (2018) (https://www.nature.com/articles/nprot.2017.115).

The package handles the different stages processing fastq files to pairs of interacting RNAs and some statistics. It *does not* handle quality issues, adapter removing etc. so the fastq files should be treated with cutadapt or equivalent before applying this package.

Installation
------------
Either ``pip install rilseq`` or ``conda install rilseq``

Simple mapping
--------------
The first stage is to map the reads to the genome. In this readme file it's assumed that the reads are paired-end but all the commands will work for single-end mapping.

Run with::

    $ map_single_fragments.py  genome.fa -g annotations.gff -1 lib1_1.fastq[.gz] lib2_1.fastq[.gz] lib3_1.fastq[.gz] -2 lib1_2.fastq[.gz] lib2_2.fastq[.gz] -d output_dir -o output_head [-r] -m max_mismatches

This script uses bwa to map the reads to the genome. You can specify two files for each library or one by using -2 or just -1. The order of the libraries in -2 should be the same as in -1 but can be shorter if some of the libraries are single-end sequenced (in this case the single-ends should be the last on -1 list). There are additional parameters that can be changed, see -h for additional help.


Chimera mapping
---------------
After we have a bam file for each library (either generated using map_single_fragments.py or any other mapper) we can look for chimeric reads.

The search is done using the ends of all the fragments in the bam file. In the case of paired-end sequencing the first 25 nucleotides (or as specified by -l argument) are taken. In the case of single-end sequencing the first 25 and last 25 of each read are taken, make sure the two regions don't overlap or you won't get any results.

All the reads in the bam files are mapped to the genome and written to the results file unless -a filename or -A are defined, in the first case the single fragments will be written to the specified file. The single fragments are used for the statistical test and will be marked as *single* so the interactions will be tested according to them but they won't be tested using Fisher's exact test. If, for some reason, a fragment was a concordant mapping as single and here the single mapping was not an option it will be omitted.

The 25 nt long sequences are screened using dust filter to remove reads with low complexity. The default threshold for the dust filter is 10, it can be changed using the --dust_thr parameter, when 0 the filter won't be used.

After the two ends are extracted they are being mapped to the genome using bwa and screened again to see if they can be on the same transcript. In order to do so we allow a relatively large number of mismatches (3 by default, set with --max_mismatches) and test if any combination of the positions each read was mapped to can result from the same transcript. We remove pairs of reads that are 1000 nt apart and map in opposite directions either as expected or in reverse order (reads which result from circular RNAs, omit this option using --keep_circular). If the -t argument is given, the pairs are tested to see if they might reside from the same transcript even if they're distance is larger than 1000 nt. This option is very useful in screening rRNAs that sometimes come from long transcripts.

After reads that might be concordant are removed, the lowest position on the chromosome of each read is collected if the read is mapped with at most 1 mismatch (set with --allowed_mismatches).

For each pair of reads the output file will contain a line with the coordinates of the first read, the second read, the read name and the word chimera/single according to the nature of the fragment. All the bam files that were given as input will be joined to the same output file, they can be further separated using the read names. Alternatively, you can run each bam file separately and cat all the reads files.

Run with::

    $ map_chimeric_fragments.py genome.fa lib1_bwa.bam lib2_bwa.bam ...

Here as above, paired-end and single-end sequencing results can be used.

Reporting Over-represented Interacting Region
---------------------------------------------
Since the sequencing results contain non-specific chimeras, another stage is needed to remove the noise of the experiment. This is done by comparing the number of reads supporting an interaction to the number of reads expected at random. Simply put we compute a 2x2 contingency table with the number of reads like this:


=============  ========  =============
 #             region 1  other regions
=============  ========  =============
region 2         K            L
other regions    M            N
=============  ========  =============

L and M include the number of single reads as well, the statistics test if two regions interact more frequently than expected by random.
If K is larger than expected the two regions are probably in actual interaction
*in vivo*. The odds ratio is computed by (K/L)/(M/N) and the p-value is computed using Fisher's exact test, testing only if K is larger than expected (single-tailed test)

The interacting regions are searched by dividing the genome to 100 nt non-overlapping windows (set with --seglen) and testing if the number of interactions is larger than expected if there are at least 5 interactions between them (set with --min_int). If the p-value is smaller than 0.05 (--max_pv) report this pair. To avoid re-reporting adjacent regions, the regions can be expanded up to 500 nts (--maxsegs) and the combination of regions with the lower p-value will be reported (after Bonferroni correction).

In some cases the interacting RNAs are found to be significant but the fraction of RNAs bound to Hfq relative to all RNAs in the cell is small, this will mean that the interaction will have no or little effect. To correct for this problem an approximation of the fraction bound by Hfq can be calculated if the results of a total RNA experiment will be given. Use --total_RNA followed by a comma separated list of indexed bam files. The script will calculate the number of reads from total RNA for each region and will estimate the abundance on Hfq, the odds ratio will be multiplied by these two values for both RNAs and reported under the "Total normalized odds ratio" column.

After the regions were determined they are decorated with additional data like genes in the region, if this is a known target, the counts of single fragments in the two genes etc.

Run this with::

     $ RILseq_significant_regions.py reads_in_list --ec_dir EcoCyc_dir --EC_chrlist "chr,COLI-K12" -t known_targets_file -c single_counts_file -r REP_elements_table

There are more arguments, some mentioned above, other can be seen using -h. In order to get gene annotations you should get the EcoCyc flat files of your organism, they require registration, point the data directory with --ec_dir. The names of the chromosomes are probably different from the bam file (the genome.fa file you used for mapping) and the EcoCyc files. You can give the script a dictionary from the bam to EcoCyc using a comma separated list of names where the name in EcoCyc follows the name in the bam file.

In addition to printing the interactions, this script can compute the interaction free-energy using RNAup (version 1 only, version 2 doesn't work) if --shuffles is > 0, it uses shuffled sequences to compute a p-value on this energy. 


Generating Plots and Tracks
---------------------------
The script plot_circos_plot read the output of map_chimeric_fragments.py to 
generate a list of interactions between regions in the chromosome. It can't 
show interactions between two chromosomes.

Together with the conf files in the data/E_coli_K12 dir and the short script
plot_interactions.sh found in this directory you can plot the interactions
with the sRNAs, rRNAs and tRNAs on the genome.

You should execute plot_interactions.sh from the directory it resides in or
give the path to the conf files. run::

    plot_interactions.sh interactions.txt interactions_plot.png

(other formats are also available like svg)

The coverage of single fragments can be viewed in UCSC genome browser for instance using the wiggle file generated by map_single_reads.py. The reads of the chimeric fragments can be written to a bed file using generate_BED_file_of_endpoints.py. The file print the position of each read in a bam file that was found to be chimeric. There is an option to print only the fragments that are part of a significant interaction, use -s interactions_file.txt to do it. When using -s you can specify a gene name (an EcoCyc ID) and generate a bed file with fragments that one of their side is mapped to the gene (-e ID). run generate_BED_file_of_endpoints.py -h for complete documentation.

Data Files
----------
This package works well for E. coli K12 (RefSeq NC_000913.2 genome and RefSeq NC_000913.3 genome). The data
directory contains two separate sub directories termed ver2 and ver3 for each of the two genome versions which
includes the genome \*.fa, the EcoCyc genes gff file and the EcoCyc transcripts gff file. These files and others in the ver2 and ver3 directories are based on EcoCyc version 19.0 and 20.0 respectively and include data from BioCyC(TM) pathway/genome database under license from SRI international. 
The genome should be indexed using bwa index genome.fa before using it. The two gff files can be generated using the scripts::

    generate_transcripts_gff.py EcoCyc_data_dir

and::

    generate_genes_gff.py EcoCyc_data_dir

There are two additional files in the ver2 data directory: a curated list of targets
taken from EcoCyc with slight changes and a table of REP elements (used for annotation of results), this table was downloaded from:  http://ecocyc.org/group?id=biocyc14-8223-3640227683 

Requirements
------------
This package requires
 - samtools (tested on version 1.2)
 - bwa (tested on version 0.7.12)
 - pysam
 - numpy & scipy
 - biopython

The project is hosted on github: https://github.com/asafpr/RILseq
