# RUN DIMPLE
# script for usage with command line

import argparse
from DIMPLE.DIMPLE import align_genevariation, print_all, post_qc, addgene, DIMPLE, generate_DMS_fragments
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="DIMPLE: Deep Indel Missense Programmable Library Engineering")
parser.add_argument('-wDir', help='Working directory for fasta files and output folder')
parser.add_argument('-geneFile', required=True, help='Input all gene sequences including backbone in a fasta format. Place all in one fasta file. Name description can include start and end points (>gene1 start:1 end:2)')
#parser.add_argument('-handle', default='AGCGGGAGACCGGGGTCTCTGAGC', help='Genetic handle for domain insertion. This is important for defining the linker. Currently uses BsaI (4 base overhang), but this can be swapped for SapI (3 base overhang).')
parser.add_argument('-matchSequences', action='store_const', const='match', default='nomatch', help='Find similar sequences between genes to avoid printing the same oligos multiple times. Default: No matching')
parser.add_argument('-oligoLen', required=True, type=int, help='Synthesized oligo length')
parser.add_argument('-fragmentLen', default=[], type=int, help='Maximum length of gene fragment')
parser.add_argument('-overlap', default=3, type=int, help='Enter number of bases to extend each fragment for overlap. This will help with insertions close to fragment boundary')
parser.add_argument('-DMS', action='store_const', const=True, default=False, help='Choose if you will run deep deep mutation scan')
parser.add_argument('-usage', default='human', help='Default is "human". Or select "ecoli"')
parser.add_argument('-insertions', default=False, nargs='+', help='Enter a list of insertions (nucleotides) to make at every position. Note, you should enter multiples of 3 nucleotides to maintain reading frame')
parser.add_argument('-deletions', default=False, nargs='+', help='Enter a list of deletions (number of nucleotides) to symmetrically delete (it will make deletions in multiples of 2x). Note you should enter multiples of 3 to maintain reading frame')
parser.add_argument('-include_substitutions', default=True, help='If you are running DMS but only want to insert or delete AA')
parser.add_argument('-barcode_start', default=0, help='To run DIMPLE multiple times, you will need to avoid using the same barcodes. This allows you to start at a different barcode.')
parser.add_argument('-restriction_sequence', default='CGTCTC', help='Recommended using BsmBI - CGTCTC or BsaI - GGTCTC')
parser.add_argument('-avoid_sequence', nargs='+', default=['CGTCTC', 'GGTCTC'], help='Avoid these sequences in the backbone - BsaI and BsmBI. For multiple sequnces use a space between inputs. Example -avoid_sequence CGTCTC GGTCTC')
args = parser.parse_args()

if args.wDir is None:
    if '/' in args.geneFile:
        args.wDir = args.geneFile.rsplit('/', 1)[0]+'/'
        args.geneFile = args.geneFile.rsplit('/', 1)[1]
    else:
        args.wDir = ''

#if any([x not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] for x in args.handle]):
#    raise ValueError('Genetic handle contains non nucleic bases')

#DIMPLE.handle = args.handle
DIMPLE.synth_len = args.oligoLen
if args.fragmentLen:
    DIMPLE.maxfrag = args.fragmentLen
else:
    DIMPLE.maxfrag = args.oligoLen - 62 - args.overlap  # 62 allows for cutsites and barcodes

if args.DMS:
    DIMPLE.dms = True
else:
    DIMPLE.dms = False

#  adjust primer primerBuffer
DIMPLE.primerBuffer += args.overlap

DIMPLE.avoid_sequence = args.avoid_sequence
DIMPLE.barcodeF = DIMPLE.barcodeF[int(args.barcode_start):]
DIMPLE.barcodeR = DIMPLE.barcodeR[int(args.barcode_start):]
DIMPLE.cutsite = Seq(args.restriction_sequence)
DIMPLE.avoid_sequence = [Seq(x) for x in args.avoid_sequence]
if DIMPLE.dms:
    if args.usage:
        DIMPLE.usage = args.usage
else:
    DIMPLE.usage = None

OLS = addgene(args.wDir+'/'+args.geneFile)

if args.matchSequences == 'match':
    align_genevariation(OLS)
if args.deletions:
    args.deletions = [int(x) for x in args.deletions]
if not any([DIMPLE.dms, args.insertions, args.deletions]):
    raise ValueError("Didn't select any mutations to generate")
generate_DMS_fragments(OLS, args.overlap, args.overlap, DIMPLE.dms, args.insertions, args.deletions, args.wDir)

post_qc(OLS)
print_all(OLS, args.wDir)
