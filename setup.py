from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
    
setup(name='RILseq',
      version='0.49',
      description='Processing RILSeq experiments results',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      scripts=[
        'bin/map_chimeric_fragments.py', 
	'bin/map_single_fragments.py',
        'bin/RILseq_significant_regions.py',
        'bin/generate_genes_gff.py',
        'bin/generate_transcripts_gff.py',
        'bin/plot_circos_plot.py',
        'bin/generate_BED_file_of_endpoints.py',
        'bin/count_chimeric_reads_per_gene.py',
        'bin/get_sequences_for_meme.py',
        'bin/plot_regions_interactions.py'],
      url='http://github.com/asafpr/RILseq',
      author='Asaf Peer',
      author_email='asafpr@gmail.com',
      license='MIT',
      packages=['RILseq'],
      install_requires=[
        'scipy', 'numpy', 'pysam', 'biopython'],
      include_package_data=True,
      zip_safe=False)
