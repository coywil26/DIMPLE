# DIMPLE: Deep Indel Missense Programmable Library Egineering
### Protein domain insertion via programmed oligo libraries
Python script for generating oligo libraries and PCR primers for programmed domain insertion

# Python Dependencies
biopython <br />
numpy

# Data Formats
Targeted genes must be in fasta format and include a minimum of 30 bases surrounding gene.
Entire plasmid sequence is advised to search for nonspecific amplification <br />
Final output is fasta format. One file for oligo pools and one file for PCR primers

# Position arguments
Gene start is defined as base number of first base in first codon and gene end is defined as base number of last base in last codon.
(Program will subtract 1 from gene start for python numbering)

To define gene position within given fasta file, add start:# end:# to fasta description. (Otherwise use command line prompts) <br />
'>geneA start:11 end:40'

# Running Test
Domain Insertion Scanning:
python3 run_spine.py -wDir tests -geneFile combined_fasta.fa -oligoLen 230 -mutationType DIS

Deep Mutational Scanning:
python3 run_spine.py -wDir tests -geneFile Kir.fa -oligoLen 230 -mutationType DMS -usage ecoli

# Usage
```
optional arguments:
-h, --help                 show this help message and exit
-wDir WDIR                 Working directory for fasta files and output folder
-geneFile GENEFILE         Input all gene sequences including backbone in a fasta
                           format. Place all in one fasta file. Name description
                           can include start and end points (>gene1 start:1
                           end:2)
-handle HANDLE             Genetic handle for domain insertion.  This is important
                           for defining the linker. Currently uses BsaI (4 base
                           overhang), but this can be swapped for SapI (3 base
                           overhang).
-matchSequences            Find similar sequences between genes to avoid printing
                           the same oligos multiple times. Default: No matching
-oligoLen OLIGOLEN         Synthesized oligo length
-fragmentLen FRAGMENTLEN   Maximum length of gene fragment
-overlap OVERLAP           Enter number of bases to extend each fragment for
                           overlap. This could help with insertion coverage close to
                           fragment boundary. Overlap does not add additional
                           insertions and thus no additional oligos.
-mutationType              Run deep insertion scan "DIS" or deep mutation scan "DMS"
-usage USAGE               Default is "human". Or select "ecoli"
```
