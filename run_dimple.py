# RUN DIMPLE
# script for usage with command line

import argparse
from DIMPLE.DIMPLE import align_genevariation, print_all, post_qc, addgene, DIMPLE, generate_DMS_fragments
from DIMPLE.utilities import parse_custom_mutations
from Bio.Seq import Seq
import os
import ast
import re
import warnings

import logging
from datetime import datetime

# Set up logging
logger = logging.getLogger(__name__)
log_file = "logs/DIMPLE-{:%Y-%m-%d-%s}.log".format(datetime.now())
# If log folder does not exist, create it
if not os.path.exists("logs"):
    os.makedirs("logs")
logger.basicConfig = logging.basicConfig(filename=log_file, level=logging.INFO)

logger.info("Started logging")

parser = argparse.ArgumentParser(description="DIMPLE: Deep Indel Missense Programmable Library Engineering")
parser.add_argument('-wDir', help='Working directory for fasta files and output folder')
parser.add_argument('-geneFile', required=True, help='Input all gene sequences including backbone in a fasta format. Place all in one fasta file. Name description can include start and end points (>gene1 start:1 end:2)')
parser.add_argument('-handle', default='AGCGGGAGACCGGGGTCTCTGAGC', help='Genetic handle for domain insertion. This is important for defining the linker. Currently uses BsaI (4 base overhang), but this can be swapped for SapI (3 base overhang).')
parser.add_argument('-dis', default=False, help='use the handle to insert domains at every position in POI')
parser.add_argument('-matchSequences', action='store_const', const='match', default='nomatch', help='Find similar sequences between genes to avoid printing the same oligos multiple times. Default: No matching')
parser.add_argument('-oligoLen', type=int, default=230, help='Synthesized oligo length')
parser.add_argument('-fragmentLen', default=0, type=int, help='Maximum length of gene fragment')
parser.add_argument('-overlap', default=4, type=int, help='Enter number of bases to extend each fragment for overlap. This will help with insertions close to fragment boundary')
parser.add_argument('-DMS', action='store_const', const=True, default=False, help='Choose if you will run deep deep mutation scan')
parser.add_argument('-custom_mutations', default=None, help='Path to file that includes custom mutations with the format position:AA')
parser.add_argument('-usage', default='human', help='Default is "human". Or select "ecoli. Or change code"')
parser.add_argument('-insertions', default=False, nargs='+', help='Enter a list of insertions (nucleotides) to make at every position. Note, you should enter multiples of 3 nucleotides to maintain reading frame')
parser.add_argument('-deletions', default=False, nargs='+', help='Enter a list of deletions (number of nucleotides) to symmetrically delete (it will make deletions in multiples of 2x). Note you should enter multiples of 3 to maintain reading frame')
parser.add_argument('-include_substitutions', default=False, help='If you are running DMS but only want to insert or delete AA')
parser.add_argument('-barcode_start', default=0, help='To run DIMPLE multiple times, you will need to avoid using the same barcodes. This allows you to start at a different barcode.')
parser.add_argument('-restriction_sequence', default='CGTCTC(G)1/5', help='Recommended using BsmBI - CGTCTC(G)1/5 or BsaI - GGTCTC(G)1/5. Do not use N')
parser.add_argument('-avoid_sequence', nargs='+', default=['CGTCTC', 'GGTCTC'], help='Avoid these sequences in the backbone - BsaI and BsmBI. For multiple sequnces use a space between inputs. Example -avoid_sequence CGTCTC GGTCTC')
parser.add_argument('-include_stop_codons', help='Include stop codons in the list of scanning mutations.', default=False, const=True, action='store_const')
parser.add_argument('-include_synonymous', help='Include synonymous codons in the list of scanning mutations.', default=False, const=True, action='store_const')
parser.add_argument('-make_double', help='Make each combination of mutations within a fragment', default=False, const=True, action='store_const')
parser.add_argument('-maximize_nucleotide_change', help='Maximize the number of nucleotide changes in each codon for easier detection in NGS', default=False, const=True, action='store_const')
parser.add_argument("-seed", help="Seed for random number generation", default=None)
args = parser.parse_args()

if args.wDir is None:
    if '/' in args.geneFile:
        args.wDir = args.geneFile.rsplit('/', 1)[0]+'/'
        args.geneFile = args.geneFile.rsplit('/', 1)[1]
    else:
        args.wDir = ''

#if any([x not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] for x in args.handle]):
#    raise ValueError('Genetic handle contains non-nucleic bases')
DIMPLE.handle = args.handle


DIMPLE.synth_len = args.oligoLen
overlapL = int(args.overlap)
overlapR = int(args.overlap)

if args.deletions:
    overlapR = max([int(x) for x in args.deletions]) + overlapR - 3

if args.fragmentLen != 0:
    DIMPLE.maxfrag = args.fragmentLen
    logger.info(f'Maximum fragment length: {DIMPLE.maxfrag} based on input')
else:
    DIMPLE.maxfrag = args.oligoLen - 64 - overlapL - overlapR # 64 allows for cutsites and barcodes
    logger.info(f'Maximum fragment length: {DIMPLE.maxfrag} based on oligo length and overlap: 2 * {args.overlap}')

#  adjust primer primerBuffer
DIMPLE.primerBuffer += overlapL

DIMPLE.avoid_sequence = args.avoid_sequence
DIMPLE.barcodeF = DIMPLE.barcodeF[int(args.barcode_start):]
DIMPLE.barcodeR = DIMPLE.barcodeR[int(args.barcode_start):]

# Check whether restriction sequence specified as enzyme or sequence
if re.match(r'[ACGT]+\([ACGT]\)\d+/\d+', args.restriction_sequence):
    tmp_cutsite = args.restriction_sequence.split('(')
    DIMPLE.cutsite = Seq(tmp_cutsite[0])
    DIMPLE.cutsite_buffer = Seq(tmp_cutsite[1].split(')')[0])
    tmp_overhang = tmp_cutsite[1].split(')')[1].split('/')
    DIMPLE.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
    if DIMPLE.cutsite == Seq('GGTCTC') and DIMPLE.cutsite_buffer == Seq('G') and DIMPLE.cutsite_overhang == 4:
        DIMPLE.enzyme = 'BsaI'
    elif DIMPLE.cutsite == Seq('CGTCTC') and DIMPLE.cutsite_buffer == Seq('G') and DIMPLE.cutsite_overhang == 4:
        DIMPLE.enzyme = 'BsmBI'
    else:
        DIMPLE.enzyme = None
elif args.restriction_sequence.upper() in ['BSAI', 'BSMBI']:
    if args.restriction_sequence.upper() == 'BSAI':
        DIMPLE.cutsite = Seq('GGTCTC')
        DIMPLE.cutsite_buffer = Seq('G')
        DIMPLE.cutsite_overhang = 4
        DIMPLE.enzyme = 'BsaI'
    else:
        DIMPLE.cutsite = Seq('CGTCTC')
        DIMPLE.cutsite_buffer = Seq('G')
        DIMPLE.cutsite_overhang = 4
        DIMPLE.enzyme = 'BsmBI'

else:
    raise ValueError(f'Restriction sequence {args.restriction_sequence} not recognized. Please check input.')

DIMPLE.avoid_sequence = [Seq(x) for x in args.avoid_sequence]

# Check whether restriction sequence is included in the avoid list
if DIMPLE.cutsite not in DIMPLE.avoid_sequence:
    DIMPLE.avoid_sequence.append(DIMPLE.cutsite)
    warnings.warn(f'Restriction sequence {DIMPLE.cutsite} was not included in the avoid list. Adding before continuing.')
    logger.warning(f'Restriction sequence {DIMPLE.cutsite} was not included in the avoid list. Adding before continuing.')

# Set up DMS parameters
DIMPLE.dms = args.DMS
DIMPLE.stop_codon = args.include_stop_codons
DIMPLE.make_double = args.make_double
DIMPLE.maximize_nucleotide_change = args.maximize_nucleotide_change

if args.custom_mutations:
    # load file with custom mutations
    with open(args.custom_mutations) as f:
        custom_mutations = f.readlines()
    # parse custom mutations
    custom_mutations = parse_custom_mutations(custom_mutations)
else:
    custom_mutations = None

# Set up random seed
if args.seed:
    DIMPLE.random_seed = int(args.seed)
else:
    DIMPLE.random_seed = None


if args.usage == 'ecoli':
    DIMPLE.usage = {
        'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.09,
        'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
        'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
        'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
        'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.3, 'TGG': 1,
        'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
        'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
        'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15
    }  # E.coli codon usage table
elif args.usage == 'human':
    DIMPLE.usage = {
        'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.2,
        'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
        'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1, 'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
        'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
        'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1,
        'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
        'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2,
        'GCT': 0.26, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25
    }
else:
    with open(args.usage) as f:
        usage = f.readlines()
    DIMPLE.usage = ast.literal_eval(usage.strip('\n'))

OLS = addgene(os.path.join(args.wDir, args.geneFile).strip())

if args.matchSequences == 'match':
    align_genevariation(OLS)
if args.deletions:
    args.deletions = [int(x) for x in args.deletions]
if not any([DIMPLE.dms, args.insertions, args.deletions]):
    raise ValueError("Didn't select any mutations to generate")

logger.info('Generating DMS fragments')

generate_DMS_fragments(OLS, overlapL, overlapR, args.include_synonymous, custom_mutations, DIMPLE.dms, args.insertions, args.deletions, args.dis, args.wDir)

post_qc(OLS)
print_all(OLS, args.wDir)

logger.info('Finished')
