# RUN SPINE
# script for usage with command line

import argparse
from SPINE.SPINE import align_genevariation, generate_DIS_fragments, print_all, post_qc, addgene, SPINEgene

parser = argparse.ArgumentParser(description="SPINE Saturated Programmable INsertion Engineering")
parser.add_argument('-wDir', help='Working directory for fasta files and output folder')
parser.add_argument('-geneFile', required=True, help='Input all gene sequences including backbone in a fasta format. Place all in one fasta file. Name description can include start and end points (>gene1 start:1 end:2)')
parser.add_argument('-handle', default='AGCGGGAGACCGGGGTCTCTGAGC', help='Genetic handle for domain insertion')
parser.add_argument('-matchSequences', action='store_const', const='match', default='nomatch', help='Find similar sequences between genes to avoid printing the same oligos multiple times. Default: No matching')
parser.add_argument('-oligoLen', required=True, type=int)
parser.add_argument('-fragmentLen', default=[], type=int)
#parser.add_argument('-restrictionSeq', default=['GGTCTC', 'CGTCTC', 'GCTCTTC'])  # BsaI, BsmBI, SapI
args = parser.parse_args()


if args.wDir is None:
    if '/' in args.geneFile:
        args.wDir = args.geneFile.rsplit('/', 1)[0]+'/'
        args.geneFile = args.geneFile.rsplit('/', 1)[1]
    else:
        args.wDir = ''

if any([x not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] for x in args.handle]):
    raise ValueError('Genetic handle contains non nucleic bases')

SPINEgene.handle = args.handle
SPINEgene.synth_len = args.oligoLen
if args.fragmentLen:
    SPINEgene.maxfrag = args.fragmentLen
#SPINEgene.restrict_seq = args.restrictionSeq

OLS = addgene(args.wDir+'/'+args.geneFile)

if args.matchSequences == 'match':
    align_genevariation(OLS)

generate_DIS_fragments(OLS, args.wDir)

post_qc(OLS)
print_all(OLS, args.wDir)
