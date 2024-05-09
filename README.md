[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/coywil26/DIMPLE/blob/master/DIMPLE.ipynb)
# DIMPLE: Deep Indel Missense Programmable Library Engineering
### Protein domain insertion via programmed oligo libraries
Python script for generating oligo libraries and PCR primers for programmed domain insertion <br />
Paper is linked here: (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02880-6)
and here: (https://www.biorxiv.org/content/10.1101/2022.07.26.501589v1.full.pdf)

# License

This code is licensed under the terms of the MIT license: 
[License](License.txt)

# Python Dependencies
biopython <br />
numpy

# Data Formats
Targeted genes must be in fasta format and include the entire plasmid sequence to search for nonspecific amplification <br />
Final output is fasta format. One file for oligo pools and one file for PCR primers

# Position arguments
Gene start is defined as base number of first base in first codon and gene end is defined as base number of last base in last codon.
(Program will subtract 1 from gene start for python numbering)

To define gene position within given fasta file, add start:# end:# to fasta description. (Otherwise use command line prompts) <br />
'>geneA start:11 end:40'

# Running Test
Deep Mutational Scanning with amino acid subsitution with stop codons and ecoli codon usage, GGG insertion, 3 base deletion, and BsaI golden gate: <br />
python run_dimple.py -wDir tests -geneFile combined_fasta.fa -oligoLen 230 -DMS -usage ecoli -include_stop_codons -restriction_sequence GGTCTC -avoid_sequence GGTCTC -insertions GGG -deletions 3

# Commandline Usage
```
arguments:
  -h, --help                                  show this help message and exit
  -wDir WDIR                                  Working directory for fasta files and output folder
  -geneFile GENEFILE                          Input all gene sequences including backbone in a fasta format. Place all in one fasta file. 
                                              Name description can include start and end points (>gene1 start:1 end:2)
  -matchSequences                             Find similar sequences between genes to avoid printing the same oligos multiple times. Default: No matching
  -oligoLen OLIGOLEN                          Synthesized oligo length
  -fragmentLen FRAGMENTLEN                    Maximum length of gene fragment
  -overlap OVERLAP                            Enter number of bases to extend each fragment for overlap. This will help with insertions close to fragment boundary
  -DMS                                        Choose if you will run deep deep mutation scan
  -usage USAGE                                Default is "human". Or select "ecoli. Or change code"
  -insertions INSERTIONS                      Enter a list of insertions (nucleotides) to make at every position. 
                                              Note, you should enter multiples of 3 nucleotides to maintain reading frame
  -deletions DELETIONS                        Enter a list of deletions (number of nucleotides) to symmetrically delete (it will make deletions in multiples of 2x). 
                                              Note you should enter multiples of 3 to maintain reading frame
  -barcode_start BARCODE_START                To run DIMPLE multiple times, you will need to avoid using the same barcodes. 
                                              This allows you to start at a different barcode.
  -restriction_sequence RESTRICTION_SEQUENCE  Recommended using BsmBI - CGTCTC or BsaI - GGTCTC
  -avoid_sequence AVOID_SEQUENCE              Avoid these sequences in the backbone - BsaI and BsmBI. For multiple sequnces use a space between inputs. 
                                              Example -avoid_sequence CGTCTC GGTCTC
  -include_stop_codons                        Include stop codons in the list of scanning mutations.
  -include_synonymous                         Include synonymous codons in the list of scanning mutations.
```

# GUI Usage
```
Run the GUI with the commandline prompt: python run_dimple_gui.py
```
![DIMPLE_GUI](https://user-images.githubusercontent.com/25623801/229153990-dc4a82e8-31f7-4914-b078-fbdb8c855761.png)
```
You must select:
  -Target Gene File - select fasta file with name description including start and end points (>gene1 start:1 end:2)
  -One or more of the mutations to make to the target gene
  
 Run the program by pressing Run DIMPLE
```
