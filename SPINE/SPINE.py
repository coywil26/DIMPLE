"""
SPINEgene Saturated Programmable INsertion Engineering

Python 3.7 package for generating oligo fragments with respective primers for saturated domain insertion for any gene of interest

Written By: David Nedrud

Requires installation of Biopython
Simple installation command: pip install biopython

File input must be .fasta/.fa format and must include the whole plasmid for primer specificity and binding
File output will also be .fasta format

Genes with variable sections can be aligned to save library space (avoid synthesizing the same sequence multiple times)
Use align_genevariation()

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio import pairwise2
import numpy as np
import os
from math import ceil
import itertools
from difflib import SequenceMatcher
# from Bio import AlignIO
# from Bio.Align.Applications import MafftCommandline

"""Synthesized Fragment Insertion"""


# This allows for individual fasta files or fasta files with multiple genes


def addgene(genefile, start=[], end=[]):
    if genefile[-1] == ' ':
        genefile = genefile[0:-1]
    tmpgene = list(SeqIO.parse(genefile.replace('\\', ''), "fasta"))
    tmpOLS = []
    for gene in tmpgene:
        if 'start:' in gene.description and 'end:' in gene.description:
            start = int(gene.description.split('start:')[1].split(' ')[0]) - 1
            end = int(gene.description.split('end:')[1].split(' ')[0])
            gene.filename = genefile.replace('\\', '')
            tmpOLS.append(SPINEgene(gene, start, end))
        else:
            gene.filename = genefile.replace('\\', '')
            tmpOLS.append(SPINEgene(gene, start, end))
    return tmpOLS  # only return the class object itself if one gene is given


class SPINEgene:
    """Synthesized Domain Insertion"""

    # Calculate and update maxfrag - Max number of nucleotides that a fragment can carry
    @property
    def synth_len(self):
        return self.__breaksites

    @synth_len.setter
    def synth_len(self, value):
        self._synth_len = value
        self.maxfrag = value - 62

    # Shared variables
    # synth_len = 230  # default
    minfrag = 24  # Picked based on smallest size for golden gate fragment efficiency
    # handle = "AGCGGGAGACCGGGGTCTCTGAGC"  # Inserted handle for domain insertion

    primerBuffer = 30  # This extends the sequence beyond the ORF for a primer. Must be greater than 30
    allhangF = []
    allhangR = []
    primerTm = (55.5, 57)  # Melting temperature limits for primers
    # Load Barcodes
    dataDirectory = os.path.abspath(os.path.dirname(__file__))
    try:
        barcodeF = list(SeqIO.parse(dataDirectory + '/data/forward_finalprimers.fasta', "fasta"))
        barcodeR = list(SeqIO.parse(dataDirectory + '/data/reverse_finalprimers.fasta', "fasta"))
    except FileNotFoundError:
        raise ValueError("Could not find barcode files. Please upload your own or place standard barcodes in the data file.")

    # @property
    # def barcodeF(self):
    #     return self._barcodeF
    # @property
    # def barcodeR(self):
    #     return self._barcodeR
    # @barcodeF.setter
    # def barcodeF(self,filename):
    #     self._barcodeF = list(SeqIO.parse(filename,"fasta"))
    # @barcodeR.setter
    # def barcodeR(self,filename):
    #     self._barcodeF = list(SeqIO.parse(filename,"fasta"))

    def __init__(self, gene, start=[], end=[]):
        #  Search for ORF
        try:
            SPINEgene.maxfrag # if SPINEgene.maxfrag doesnt exist, create it
        except AttributeError:
            SPINEgene.maxfrag = self.synth_len - 62  # based on space for barcodes, cut sites, handle
        self.geneid = gene.name
        self.linked = set()
        self.genePrimer = []
        self.oligos = []
        self.barPrimer = []
        self.fullGene = gene.seq
        self.split = 0
        self.num_frag_per_oligo = 1
        self.doublefrag = 0
        self.filename = gene.filename
        # Set up variables. Could have this as user input in the class
        if self.usage == 'ecoli':
            self.usage = {
                'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.09,
                'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
                'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
                'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
                'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.3, 'TGG': 1,
                'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
                'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
                'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15
                }  # E.coli codon usage table
        if self.usage == 'human':
            self.usage = {
                'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.2,
                'CTT': 0.13, 'CTC': 0.2, 'CTA': 0.07, 'CTG': 0.41, 'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
                'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1, 'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
                'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
                'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1,
                'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11, 'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
                'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.2, 'AGG': 0.2,
                'GCT': 0.26, 'GCC': 0.4, 'GCA': 0.23, 'GCG': 0.11, 'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25
                }
        self.SynonymousCodons = {
            'Cys': ['TGT', 'TGC'],
            'Asp': ['GAT', 'GAC'],
            'Ser': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
            'Gln': ['CAA', 'CAG'],
            'Met': ['ATG'],
            'Asn': ['AAC', 'AAT'],
            'Pro': ['CCT', 'CCG', 'CCA', 'CCC'],
            'Lys': ['AAG', 'AAA'],
            'STOP': ['TAG', 'TGA', 'TAA'],
            'Thr': ['ACC', 'ACA', 'ACG', 'ACT'],
            'Phe': ['TTT', 'TTC'],
            'Ala': ['GCA', 'GCC', 'GCG', 'GCT'],
            'Gly': ['GGT', 'GGG', 'GGA', 'GGC'],
            'Ile': ['ATC', 'ATA', 'ATT'],
            'Leu': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
            'His': ['CAT', 'CAC'],
            'Arg': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
            'Trp': ['TGG'],
            'Val': ['GTA', 'GTC', 'GTG', 'GTT'],
            'Glu': ['GAG', 'GAA'],
            'Tyr': ['TAT', 'TAC']}
        self.aminoacids = ['Cys', 'Asp', 'Ser', 'Gln', 'Met', 'Asn', 'Pro', 'Lys', 'Thr', 'Phe', 'Ala', 'Gly', 'Ile', 'Leu',
                      'His', 'Arg', 'Trp', 'Val', 'Glu', 'Tyr']

        # First check for BsaI sites and BsmBI sites
        if any([gene.seq.upper().count(cut) for cut in ['GGTCTC', 'GAGACC', 'CGTCTC', 'GAGACG']]):
            raise ValueError('Unwanted Restriction cut sites found. Please input plasmids with these removed.')  # change codon
        if start and end and (end - start) % 3 != 0:
            print('Gene length is not divisible by 3')
            start = []
            end = []
        if not start and not end:
            # Scan through all strands and frames for open reading frames
            min_protein_len = 100  # Support for finding gene position in vector
            genestart = []
            geneend = []
            genestrand = []
            print("Analyzing Gene:" + self.geneid)
            for strand, nuc in [(+1, gene.seq), (-1, gene.seq.reverse_complement())]:
                for frame in range(3):
                    length = 3 * ((len(gene) - frame) // 3)  # Multiple of three
                    translated = nuc[frame:frame + length].translate()
                    for protein in translated.split("*"):
                        if len(protein) >= min_protein_len:
                            if len(protein.split('M', 1)) > 1:
                                ORF = 'M' + protein.split('M', 1)[1]
                                genestart.append(translated.find(ORF) * 3 + frame + 1)
                                geneend.append(genestart[-1] + len(ORF) * 3 - 1)
                                genestrand.append(strand)
                                print("ORF#%i %s...%s - length %i, strand %i, frame %i" % (len(genestart), ORF[:25], ORF[-10:], len(ORF), strand, frame + 1))
            # Select Gene Frame
            while True:
                try:
                    genenum = int(input("Which ORF are you targeting? (number):"))  # userinput
                    if genestrand[genenum - 1] == -1:
                        gene.seq = gene.seq.reverse_complement()
                except:
                    print("Please enter number")
                    continue
                else:
                    break
            start = genestart[genenum - 1] - 1  # subtract 1 to account for 0 indexing
            end = geneend[genenum - 1]
            print(gene.seq[start:end].translate()[:10])
            quest = "g" # holding place

            while quest != "n" and quest != "y":
                quest = input("Is this the beginning of your gene?(position %i) (y/n):" % (start))
            while quest == "n":
                start = int(input("Enter the starting position your gene:"))
                print(gene.seq[start:end].translate()[:10])
                quest = input("Is this the beginning of your gene?(position %i) (y/n):" % (start))
            print(gene.seq[start:end].translate()[-10:])
            quest = "g"
            while quest != "n" and quest != "y":
                quest = input("Is the size of your gene %ibp? (y/n):" % (end - start))
            while quest == "n":
                try:
                    end = int(input("Enter nucleotide length of your gene:")) + start
                    if (end - start) % 3 != 0:
                        print('Length is not divisible by 3')
                        continue
                    print(gene.seq[start:end].translate()[-10:])
                    quest = input("Is this end correct? (y/n):")
                except:
                    print("Please enter a number")
                    quest = "n"
        self.aacount = int((end - start) / 3)
        self.start = start
        self.end = end
        # record sequence with extra bp to account for primer. for plasmids (circular) we can rearrange linear sequence)
        if start - self.primerBuffer < 0:
            self.seq = gene.seq[start - self.primerBuffer:] + gene.seq[:end + self.primerBuffer]
        elif end + self.primerBuffer > len(gene.seq):
            self.seq = gene.seq[start - self.primerBuffer:] + gene.seq[:end + self.primerBuffer - len(gene.seq)]
        else:
            self.seq = gene.seq[start - self.primerBuffer:end + self.primerBuffer]
        # Determine Fragment Size and store beginning and end of each fragment
        num = int(round(((end - start) / float(SPINEgene.maxfrag)) + 0.499999999))  # total bins needed (rounded up)
        insertionsites = range(start + 3, end + 3, 3)
        fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]
        # if any(x<144 for x in fragsize):
        #     raise ValueError('Fragment size too low')
        print('Initial Fragment Sizes for:' + self.geneid)
        print(fragsize)
        total = SPINEgene.primerBuffer
        breaksites = [SPINEgene.primerBuffer]
        for x in fragsize:
            total += x
            breaksites.extend([total])
        self.breaklist = [[x + 3, x + fragsize[idx]] for idx, x in
                          enumerate(breaksites[:-1])]  # insertion site to insertion site
        self.problemsites = set()
        self.unique_Frag = [True]*len(fragsize)
        self.fragsize = fragsize
        self.__breaksites = breaksites

    def ochre(self):
        if len(self.SynonymousCodons['STOP']) < 2:
            raise Exception('You have removed all stop codons')
        self.usage_ecoli['TAG'] = 1
        self.usage_human['TAG'] = 1
        # SynonymousCodons['STOP'] = ['TGA','TAA']
        del self.SynonymousCodons['STOP'][0]
        self.SynonymousCodons['OCHRE'] = ['TAG']
        self.aminoacids.extend('OCHRE')

    def amber(self):
        if len(self.SynonymousCodons['STOP']) < 2:
            raise Exception('You have removed all stop codons')
        self.usage_ecoli['TAA'] = 1
        self.usage_human['TAA'] = 1
        del self.SynonymousCodons['STOP'][2]
        self.SynonymousCodons['AMBER'] = ['TAA']
        self.aminoacids.extend('AMBER')

    def opal(self):
        if len(self.SynonymousCodons['STOP']) < 2:
            raise Exception('You have removed all stop codons')
        self.usage_ecoli['TGA'] = 1
        self.usage_human['TGA'] = 1
        del self.SynonymousCodons['STOP'][1]
        self.SynonymousCodons['OPAL'] = ['TGA']
        self.aminoacids.extend('OPAL')

    def __getitem__(self):
        return

    # Update Breaksites
    @property
    def breaksites(self):
        return self.__breaksites

    @breaksites.setter
    def breaksites(self, value):
        if isinstance(value, list):
            if any([(x - SPINEgene.primerBuffer) % 3 != 0 for x in value]):
                raise ValueError('New Breaksites are not divisible by 3')
            if value[0] != self.breaksites[0] or value[-1] != self.breaksites[-1]:
                if input('Beginning and End of gene have changed. Are you sure you want to continue? (y/n)') != 'y':
                    raise Exception('Canceled user set break sites')
            self.__breaksites = value
            fragsize = [j - i for i, j in zip(value[:-1], value[1:])]
            self.fragsize = fragsize
            self.breaklist = [[x + 3, x + fragsize[idx]] for idx, x in
                              enumerate(value[:-1])]  # insertion site to insertion site
            print('New Fragment Sizes for:' + self.geneid)
            print(fragsize)
            # fragment_genes(self)

        else:
            raise ValueError('Breaklist input is not a list')


def align_genevariation(OLS):
    if not isinstance(OLS[0], SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    match = []
    print('------------Finding homologous regions------------')
    # First find genes with matching sequences
    for m in range(len(OLS)):
        remlist = range(len(OLS))[m + 1:]
        for p in remlist:
            # score = pairwise2.align.globalmx(OLS[m].seq[genedata[m][2][0][0]:genedata[m][2][-1][0]],OLS[p].seq[genedata[p][2][0][0]:genedata[p][2][-1][0]],2,-1,one_alignment_only=True)
            score = pairwise2.align.globalmx(OLS[m].seq, OLS[p].seq, 2, -1, one_alignment_only=True)
            if score[0][2] / score[0][4] > 1.5:  # Threshold for a matched gene set
                index = [x for x, geneset in enumerate(match) if m in geneset or p in geneset]  # Determine if aligned genes are in any of the previously match sets
                if not index:  # Create a new set if not
                    print(OLS[m].geneid)
                    print(OLS[p].geneid)
                    if input("Are these genes linked? (y/n):") == 'y':
                        match.append(set([m, p]))
                else:
                    if p not in match[index[0]] or m not in match[index[0]]:
                        for items in match[index[0]].union(set([p, m])):
                            print(OLS[items].geneid)
                        if input("Are these genes linked? (y/n):") == 'y':
                            match[index[0]].add(p)
                            match[index[0]].add(m)
    # Create fragments for each match
    if match:
        for tmpset in match:
            matchset = list(tmpset)
            print('Determining Gene Variation for genes:' + ','.join([OLS[i].geneid for i in matchset]))
            # Create Alignments for each gene set
            # write to file
            # with open('tmp_align_seq.fasta', 'w') as f:
            #     for geneidx in matchset:
            #         f.write('>%s %s\n' % (OLS[geneidx].geneid,str(geneidx)))
            #         f.write('%s\n' % OLS[geneidx].seq)
            # mafft_cline = MafftCommandline(input='tmp_align_seq.fasta')
            # stdout, stderr = mafft_cline()
            # with open('aligned.fasta','w') as f:
            #     f.write(stdout)
            # align = AlignIO.read('aligned.fasta','fasta')
            max_gene_len = 0
            variablesites = set()
            for i, j in itertools.combinations(matchset, 2):
                max_gene_len = max(max_gene_len, len(OLS[i].seq)-2 * SPINEgene.primerBuffer, len(OLS[j].seq)-2 * SPINEgene.primerBuffer)
                seq_match = SequenceMatcher(None, OLS[i].seq, OLS[j].seq)
                # Determine variable regions
                variablesites.update([x.size for x in seq_match.get_matching_blocks() if x.size != len(OLS[i].seq) and x.size != len(OLS[j].seq) and x.size != 0])  # not sure how to account for zero
            problemsites = set()
            for kk in variablesites:
                problemsites.update(range(kk - SPINEgene.primerBuffer, kk + SPINEgene.primerBuffer))  # Add space for primers to bind
            # Determine Fragment Size while avoiding variable regions - must be same for all genes
            num = int(round(((max_gene_len) / float(SPINEgene.maxfrag)) + 0.499999999))  # total bins needed (rounded up)
            insertionsites = range(SPINEgene.primerBuffer + 3, max_gene_len + SPINEgene.primerBuffer + 3, 3)  # all genes start with a buffer
            fragsize = [len(insertionsites[i::num]) * 3 for i in list(range(num))]
            total = SPINEgene.primerBuffer
            breaksites = [SPINEgene.primerBuffer]  # first site is always the max primer length (adjusted at beginning)
            for x in fragsize:
                total += x
                breaksites.extend([total])
            available_sites = [xsite for xsite in range(0, max_gene_len+SPINEgene.primerBuffer+1, 3) if xsite not in problemsites]
            breaksites = [site if site in available_sites else min(available_sites, key=lambda x:abs(x-site)) for site in breaksites]  # remove problemsites?
            if any(x < SPINEgene.minfrag or x > SPINEgene.maxfrag for x in fragsize):
                print(fragsize)
                raise ValueError('Fragment size too low')  # this was decided by author. could be changed
            fragsize = [j-i for i, j in zip(breaksites[:-1], breaksites[1:])]
            breaklist = [[x + 3, x + fragsize[idx]] for idx, x in
                         enumerate(breaksites[:-1])]  # insertion site to insertion site
            unique_Frag = [[] for x in range(max(matchset)+1)]  # a list of fragments that do not match
            for x in breaklist:
                sequences = [str(OLS[i].seq[x[0]:x[1]]) for i in matchset]
                index = [matchset[i] for i, x in enumerate(sequences) if i == sequences.index(x)]
                for idx in matchset:
                    if idx in index:
                        unique_Frag[idx].extend([True])
                    else:
                        unique_Frag[idx].extend([False])

            print('Finished Alignment. Fragment Sizes for combined genes:')
            print(fragsize)
            for idx in matchset:  # setting these to the same variable should link them for processing later
                OLS[idx].problemsites = problemsites  # add gap range to problemsites variable to avoid breaking in a gap
                OLS[idx].breaklist = breaklist
                OLS[idx].fragsize = fragsize
                OLS[idx].breaksites = breaksites
                OLS[idx].linked.update(matchset)
                OLS[idx].unique_Frag = unique_Frag[idx]
    else:
        print(
            'No redundant sequences found. Matching sequences may be too short or not aligned to reduce number of oligos synthesized')


def find_geneprimer(genefrag, start, end):
    # start is variable to adjust melting temperature
    # end is fixed with restriction site added
    primer = genefrag[start:end].complement() + "CTCTGCATA" # added ATA for cleavage close to end of DNA fragment
    # Check melting temperature
    # find complementary sequences
    comp = 0  # compensate for bases that align with bsmbi
    while primer.complement()[end - start + comp] == genefrag[end + comp]:
        comp += 1
    # comp += 1 # This is important for single basepair overhang
    tm2 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN2)
    tm4 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN4)
    count = 0
    while tm2 < SPINEgene.primerTm[0] or tm2 > SPINEgene.primerTm[1] or tm4 < SPINEgene.primerTm[0] or tm4 > SPINEgene.primerTm[1]:
        if tm2 < SPINEgene.primerTm[0] or tm4 < SPINEgene.primerTm[0]:
            start += -1
            primer = genefrag[start:end].complement() + "CTCTGCATA"  # BsmbI Site addition
            tm2 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN4)
        if count > 12 or start == 0:  # stop if caught in inf loop or if linker is at max (31 + 7 = 38 bases)
            break
        if tm2 > SPINEgene.primerTm[1] and tm4 > SPINEgene.primerTm[1]:
            start += 1
            primer = genefrag[start:end].complement() + "CTCTGCATA"
            # tm = mt.Tm_NN(primer[0:e-s+comp],c_seq=genefrag[s:e+comp],nn_table=mt.DNA_NN2)
            tm2 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer[0:end - start + comp], nn_table=mt.DNA_NN4)
        count += 1
    # optional - force first nucleotide to a C or G
    # while primer[0]=="T" or primer[0]=="A" or primer[0]=="t" or primer[0]=="a":
    #     s += -1
    #     primer = genefrag[s:e].complement()+"CTCTGCA"
    #     tm = mt.Tm_NN(primer[0:e-s+comp],nn_table=mt.DNA_NN2)
    # return final primer with tm
    return primer.complement().reverse_complement(), round(tm2, 1), start


def find_fragment_primer(fragment, stop):
    start = 0  # starts at maximum length
    if stop > 25:
        end = 25
    else:
        end = stop
    count = 0
    primer = fragment[start:end]
    tm2 = mt.Tm_NN(primer, nn_table=mt.DNA_NN2)  # Two methods of finding melting temperature seems more consistent
    tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
    while tm2 < SPINEgene.primerTm[0] or tm2 > SPINEgene.primerTm[1] or tm4 < SPINEgene.primerTm[0] or tm4 > SPINEgene.primerTm[1] or len(primer) < 16:
        count += 1
        if count > 12 or end > stop:  # stop if caught in inf loop or if primer is larger than the barcode
            end = stop
            primer = fragment[start:end]
            break

        if tm2 < SPINEgene.primerTm[0] or tm4 < SPINEgene.primerTm[0]:
            if start == 0:
                break
            end += 1
            primer = fragment[start:end]
            tm2 = mt.Tm_NN(primer, nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)

        if tm2 > SPINEgene.primerTm[1] or tm4 > SPINEgene.primerTm[1]:
            end += -1
            primer = fragment[start:end]
            tm2 = mt.Tm_NN(primer, nn_table=mt.DNA_NN2)
            tm4 = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
    return primer, round(tm2, 1)


def check_nonspecific(primer, fragment, point):
    non = []
    # Forward
    for i in range(len(fragment) - len(primer)):  # Scan each position
        match = [primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))]
        first = 10
        for k in range(len(match) - 3):
            if (match[k] and match[k + 1] and match[k + 3]) or (match[k] and match[k + 1] and match[k + 2]):
                first = k
                break
        if sum(match[first:]) > len(primer[first:]) * 0.8 and sum(match[first:]) > 6 and match[-1] and point != i:  # string compare - sum of matched nt is greater than 80%
            try:
                melt = mt.Tm_NN(primer[first:], c_seq=fragment[i + first:i + len(primer)].complement(),
                                nn_table=mt.DNA_NN2, de_table=mt.DNA_DE1, imm_table=mt.DNA_IMM1)
                if melt > 20:
                    print('Found non-specfic match at ' + str(i + 1) + 'bp:')
                    print(' match:' + fragment[i:i + len(primer)])
                    print('primer:' + primer + ' Tm:' + str(round(melt, 1)))
                if melt > 35:
                    non.append(True)
            except ValueError as valerr:
                print(str(valerr) + ". Please check position manually:" + str(i + 1) + " forward")
                print('Primer:' + primer)
                print('Match: ' + fragment[i:i + len(primer)])
                non.append(False)
    # Reverse
    fragment = fragment.reverse_complement()
    for i in range(len(fragment) - len(primer)):
        match = [primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))]
        first = 10
        for k in range(0, len(match) - 3, 1):
            if match[k] and match[k + 1] and match[k + 3]:
                first = k
                break
        if sum(match[first:]) > len(primer[first:]) * 0.8 and sum(match[first:]) > 6 and match[-1] and point != -i:  # string compare - sum of matched nt is greater than 80%
            try:
                melt = mt.Tm_NN(primer[first:], c_seq=fragment[i + first:i + len(primer)].complement(),
                                nn_table=mt.DNA_NN2, de_table=mt.DNA_DE1, imm_table=mt.DNA_IMM1)
                if melt > 20:
                    print('Found non-specfic match at ' + str(i + 1) + 'bp:')
                    print(' match:' + fragment[i:i + len(primer)])
                    print('primer:' + primer + ' Tm:' + str(melt))
                if melt > 35:
                    non.append(True)
            except ValueError as valerr:
                print(str(valerr) + ". Please check position manually:" + str(i + 1) + " reverse")
                print('Primer:' + primer)
                print('Match: ' + fragment[i:i + len(primer)])
                non.append(False)
    return sum(non)


def switch_fragmentsize(gene, detectedsite, OLS):
    if not isinstance(gene, SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    start = gene.start
    skip = False
    count = 0
    print('Non specific Fragment:' + str(detectedsite))
    gene.problemsites.add(gene.breaksites[detectedsite])
    if all(item == SPINEgene.maxfrag for item in gene.fragsize) or any(
            item > SPINEgene.maxfrag for item in gene.fragsize):
        tmpmax = SPINEgene.maxfrag + 3
    else:
        tmpmax = SPINEgene.maxfrag
    while True:
        if count > len(gene.breaksites):
            tmpmax = SPINEgene.maxfrag + 3
        count += 1
        # Find connecting Fragments
        if detectedsite == 0 or detectedsite == len(gene.fragsize):
            print('Issue with primer on end of gene')
            skip = True
            break
        elif gene.fragsize[detectedsite] == gene.fragsize[detectedsite - 1] and gene.fragsize[detectedsite] >= tmpmax:
            if all(item >= tmpmax for item in gene.fragsize[detectedsite + 1:]) and not all(item >= tmpmax for item in gene.fragsize[:detectedsite - 1]):
                shift = 3
                while gene.breaksites[detectedsite] + shift in gene.problemsites:
                    shift += 3
            if all(item >= tmpmax for item in gene.fragsize[:detectedsite - 1]) and not all(item >= tmpmax for item in gene.fragsize[detectedsite + 1:]):
                shift = -3
                while gene.breaksites[detectedsite] + shift in gene.problemsites:
                    shift += -3
            else:
                if detectedsite < len(gene.fragsize) / 2:  # should be based on problemsites not where it is located in the gene
                    shift = 3
                    while gene.breaksites[detectedsite] + shift in gene.problemsites:
                        shift += 3
                else:
                    shift = -3
                    while gene.breaksites[detectedsite] + shift in gene.problemsites:
                        shift += -3
        elif gene.fragsize[detectedsite] > gene.fragsize[detectedsite - 1]:
            shift = 3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift += 3
        elif gene.fragsize[detectedsite] < gene.fragsize[detectedsite - 1]:
            shift = -3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift += -3
        elif gene.fragsize[detectedsite] == gene.fragsize[detectedsite - 1] and gene.fragsize[detectedsite] < tmpmax:
            shift = -3
            while gene.breaksites[detectedsite] + shift in gene.problemsites:
                shift = -shift
                if shift < 0:
                    shift += -3
        # Process shift and reprocess fragments
        gene.breaksites[detectedsite] = gene.breaksites[detectedsite] + shift
        gene.fragsize = [j - i for i, j in zip(gene.breaksites[:-1], gene.breaksites[1:])]
        gene.breaklist = [[x + 3, x + gene.fragsize[idx]] for idx, x in enumerate(gene.breaksites[:-1])]
        # recheck for size limit issues
        tmpsite = [topidx for topidx, item in enumerate(gene.fragsize) if item > tmpmax]
        if tmpsite:
            # pick which side to adjust
            if tmpsite[0] == len(gene.fragsize):
                detectedsite = tmpsite[0]
            elif tmpsite[0] == 0:
                detectedsite = tmpsite[0] + 1
            elif tmpsite[0] == detectedsite and tmpsite[0] + 1 < len(gene.fragsize):
                detectedsite = tmpsite[0] + 1
            else:
                detectedsite = tmpsite[0]
        else:
            break
    print(gene.fragsize)
    # align all linked genes to the same breaksites
    for tmp in gene.linked:
        OLS[tmp].breaksites = gene.breaksites
        OLS[tmp].fragsize = gene.fragsize
        OLS[tmp].breaklist = gene.breaklist
    return skip


def check_overhangs(gene, OLS, overlap):
    # Force all overhangs to be different within a gene (no more than 2 matching in a row)
    if not isinstance(gene, SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    while True:
        overhang = []
        for idx, y in enumerate(gene.breaklist):
            overhang.append([gene.seq[y[0] - 4 - overlap : y[0] - overlap], idx])  # Forward overhang
            overhang.append([gene.seq[y[1] + overlap : y[1] + 4 + overlap], idx + 1])  # Reverse overhang
        detectedsites = set()  # stores matching overhangs
        for i in range(len(overhang)):  # check each overhang for matches
            for j in [x for x in range(len(overhang)) if x != i]:  # permutate over every overhang combination to find matches
                if overhang[i][0] == overhang[j][0] or \
                        overhang[i][0][:3] == overhang[j][0][:3] or \
                        overhang[i][0][1:] == overhang[j][0][1:] or \
                        overhang[i][0] == overhang[i][0].reverse_complement():  # no matching overhangs or 3 base match or palindromes
                    detectedsites.update([overhang[i][1]])
        for detectedsite in detectedsites:
            if detectedsite == 0:
                detectedsite = 1
            print("------------------ Fragment size swapped due to matching overhangs ------------------")
            skip = switch_fragmentsize(gene, detectedsite, OLS)
        else:
            break


def generate_DIS_fragments(OLS, overlap, folder=''):
    if not isinstance(OLS[0], SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    # Loop through each gene or gene variation
    finishedGenes = []
    for ii, gene in enumerate(OLS):
        print('--------------------------------- Analyzing Gene:' + gene.geneid + ' ---------------------------------')
        if not any([tmp in finishedGenes for tmp in gene.linked]):  # only run analysis for one of the linked genes
            # Quality Control for overhangs from the same gene
            check_overhangs(gene, OLS, overlap)
        # Generate oligos and Primers
        idx = 0 # index for fragment
        totalcount = 0
        #storage for unused barcodes
        compileF = []
        compileR = []

        while idx < len(gene.breaklist):
            if idx == 0:
                gene.oligos = []
                gene.barPrimer = []
                gene.genePrimer = []
            frag = gene.breaklist[idx]
            fragstart = str(int((frag[0] - SPINEgene.primerBuffer) / 3))
            fragend = str(int((frag[1] - SPINEgene.primerBuffer) / 3))
            print('Creating Gene:' + gene.geneid + ' --- Fragment:' + fragstart + '-' + fragend)
            if not any([tmp in finishedGenes for tmp in gene.linked]):  # only run analysis for one of the linked genes
                # Primers for gene amplification with addition of BsmBI site
                genefrag = gene.seq[frag[0]-SPINEgene.primerBuffer : frag[0]+SPINEgene.primerBuffer]
                reverse, tmR, sR = find_geneprimer(genefrag, 15, SPINEgene.primerBuffer+1-overlap) # 15 is just a starting point
                genefrag = gene.seq[frag[1]-SPINEgene.primerBuffer: frag[1]+SPINEgene.primerBuffer]
                forward, tmF, sF = find_geneprimer(genefrag.reverse_complement(), 15, SPINEgene.primerBuffer+1-overlap)
                tmpr = check_nonspecific(reverse, gene.seq, frag[0] - len(gene.seq) + 10 - overlap) # negative numbers look for reverse primers
                tmpf = check_nonspecific(forward, gene.seq, frag[1] - 10 + overlap)
                if tmpf or tmpr:
                    # swap size with another fragment
                    if tmpf:
                        idx = idx + 1
                    # swap size with another fragment
                    print("------------------ Fragment size swapped due to non-specific primers ------------------")
                    skip = switch_fragmentsize(gene, idx, OLS)
                    if skip:
                        continue
                    # Quality Control for overhangs from the same gene
                    # check_overhangs(gene, OLS)
                    SPINEgene.barcodeF.extend(compileF)  # return unused primers
                    SPINEgene.barcodeR.extend(compileR)
                    compileF = []  # reset unused primers
                    compileR = []
                    gene.genePrimer = []  # reset gene all primers due to nonspecific primer
                    gene.barPrimer = []
                    idx = 0
                    continue
                # Store
                gene.genePrimer.append(SeqRecord(reverse, id=gene.geneid + "_geneP_DI-" + str(idx + 1) + "_R",
                                                 description="Frag" + fragstart + "-" + fragend + ' ' + str(tmR) + 'C'))
                gene.genePrimer.append(SeqRecord(forward, id=gene.geneid + "_geneP_DI-" + str(idx + 1) + "_F",
                                                 description="Frag" + fragstart + "-" + fragend + ' ' + str(tmF) + 'C'))
            if gene.unique_Frag[idx]:
                # Create gene fragments with insertions
                tmF = 0
                tmR = 0
                count = 0
                tmpseq = gene.seq[frag[0] - 4-overlap : frag[1]+4+overlap].ungap('-') #4 is overhang for BsmBI
                while tmF < SPINEgene.primerTm[0] or tmR < SPINEgene.primerTm[0]:  # swap out barcode if tm is low
                    difference = (SPINEgene.synth_len - (len(tmpseq) + 14 + len(SPINEgene.handle)))  # 14 bases is the length of the restriction sites with overhangs (7 bases each)
                    barF = SPINEgene.barcodeF.pop(0)
                    barR = SPINEgene.barcodeR.pop(0)
                    count += 1  # How many barcodes used
                    compileF.append(barF)
                    compileR.append(barR)
                    while difference / 2 > len(barF):
                        tmpF = SPINEgene.barcodeF.pop(0)
                        tmpR = SPINEgene.barcodeR.pop(0)
                        compileF.append(tmpF)
                        compileR.append(tmpR)
                        barF += tmpF
                        barR += tmpR
                        count += 1  # How many barcodes used
                    # check barcode for BsmBI cut sites (not BsaI)
                    if barF.seq[0:int(difference / 2)].count('CGTCTC') or barF.seq[0:int(difference / 2)].count('GAGACG') or barR.seq.reverse_complement()[0:difference - int(difference / 2)].count('CGTCTC') or barR.seq.reverse_complement()[0:difference - int(difference / 2)].count('GAGACG'):
                        print('Additional restriction sites found in barcode. Replacing barcodes')
                        # return barcodes?
                        continue
                    tmpfrag = barF.seq[0:int(difference / 2)] + "CGTCTCC" + tmpseq + "GGAGACG" + barR.seq.reverse_complement()[0:difference - int(difference / 2)]
                    # primers for amplifying subpools
                    offset = int(difference / 2) + 11  # add 11 bases for type 2 restriction
                    primerF, tmF = find_fragment_primer(tmpfrag, offset)
                    primerR, tmR = find_fragment_primer(tmpfrag.reverse_complement(), (difference - offset + 22))

                # Create a gene fragment for each insertion point. Important! This loop does all the work. Change this loop to edit mutation type
                for i in range(offset + overlap, offset + gene.fragsize[idx] + overlap, 3):  # Could replace fragsize with frag[1]-frag[0]
                    xfrag = tmpfrag[0:i] + SPINEgene.handle + tmpfrag[i:]
                    if len(xfrag) < SPINEgene.synth_len:
                        print('Fragment is ' + str(len(xfrag)))
                        raise ValueError('Fragment is less than specifified oligo length')
                    # Check each cassette for more than 2 BsmBI and 2 BsaI sites
                    if (xfrag.upper()[offset:len(tmpseq)+offset].count('GGTCTC') + xfrag.upper()[offset:len(tmpseq)+offset].count('GAGACC')) > 2:
                        raise ValueError('BsaI site found within insertion fragment')
                    if (xfrag.upper()[offset:len(tmpseq)+offset].count('CGTCTC') + xfrag.upper()[offset:len(tmpseq)+offset].count('GAGACG')) > 2:
                        raise ValueError('BsmBI site found within insertion fragment')
                    gene.oligos.append(SeqRecord(xfrag,
                                                    id=gene.geneid + "_DI_" + fragstart + "-" + fragend + "_insertion" + str(int((frag[0] + i - offset - SPINEgene.primerBuffer - overlap) / 3)),
                                                    description=''))
                # Store primers for gene fragment
                gene.barPrimer.append(SeqRecord(primerF, id=gene.geneid + "_oligoP_DI-" + str(idx + 1) + "_F",
                                                description="Frag" + fragstart + "-" + fragend + "_" + str(tmF) + 'C'))
                gene.barPrimer.append(SeqRecord(primerR, id=gene.geneid + "_oligoP_DI-" + str(idx + 1) + "_R",
                                                description="Frag" + fragstart + "-" + fragend + "_" + str(tmR) + 'C'))
                print('Barcodes used:' + str(count))
                print('Barcodes Remaining:' + str(len(SPINEgene.barcodeF)))
            idx += 1
        # Export files (fasta)
        # Fragments
        SeqIO.write(gene.oligos, os.path.join(folder.replace('\\', ''),
                                                 gene.geneid + "_DI_Oligos.fasta"), "fasta")
        # Barcode Primers
        SeqIO.write(gene.barPrimer, os.path.join(folder.replace('\\', ''),
                                                 gene.geneid + "_DI_Oligo_Primers.fasta"), "fasta")
        # Amplification Primers
        SeqIO.write(gene.genePrimer, os.path.join(folder.replace('\\', ''),
                                                  gene.geneid + "_DI_Gene_Primers.fasta"), "fasta")
        # Record finished gene for aligned genes
        finishedGenes.extend([ii])


def generate_DMS_fragments(OLS, overlap, folder=''):
    if not isinstance(OLS[0], SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    # Loop through each gene or gene variation
    finishedGenes = []
    for ii, gene in enumerate(OLS):
        print('--------------------------------- Analyzing Gene:' + gene.geneid + ' ---------------------------------')
        gene.breaklist[0][0] += 3  # Do not mutate first codon
        gene.fragsize[0] += -3  # Readjust size to match breaklist
        tmpbreaklist = [[x[0] - 3, x[1]] for x in gene.breaklist]  # Shift gene fragments to mutate first codon
        if not any([tmp in finishedGenes for tmp in gene.linked]):  # only run analysis for one of the linked genes
            # Quality Control for overhangs from the same gene
            check_overhangs(gene, OLS, overlap)
        # Generate oligos and Primers
        idx = 0 # index for fragment
        totalcount = 0
        #storage for unused barcodes
        compileF = []
        compileR = []
        missingSites = []
        offset_list = []
        # missingTable = [[1]*gene.aacount]*gene.aacount
        missingFragments = []
        all_grouped_oligos = []
        while idx < len(tmpbreaklist):
            if idx == 0:
                gene.oligos = []
                gene.barPrimer = []
                gene.genePrimer = []
            frag = tmpbreaklist[idx]
            grouped_oligos = []
            fragstart = str(int((frag[0] - SPINEgene.primerBuffer) / 3) + 1)
            fragend = str(int((frag[1] - SPINEgene.primerBuffer) / 3))
            print('Creating Gene:' + gene.geneid + ' --- Fragment:' + fragstart + '-' + fragend)
            if not any([tmp in finishedGenes for tmp in gene.linked]):  # only run analysis for one of the linked genes
                # Primers for gene amplification with addition of BsmBI site
                genefrag = gene.seq[frag[0]-SPINEgene.primerBuffer : frag[0]+SPINEgene.primerBuffer]
                reverse, tmR, sR = find_geneprimer(genefrag, 15, SPINEgene.primerBuffer+1-overlap) # 15 is just a starting point
                genefrag = gene.seq[frag[1]-SPINEgene.primerBuffer: frag[1]+SPINEgene.primerBuffer]
                forward, tmF, sF = find_geneprimer(genefrag.reverse_complement(), 15, SPINEgene.primerBuffer+1-overlap)
                tmpr = check_nonspecific(reverse, gene.seq, frag[0] - len(gene.seq) + 10 - overlap) # negative numbers look for reverse primers
                tmpf = check_nonspecific(forward, gene.seq, frag[1] - 10 + overlap)
                if tmpf or tmpr:
                    # swap size with another fragment
                    if tmpf:
                        idx = idx + 1
                    # swap size with another fragment
                    print("------------------ Fragment size swapped due to non-specific primers ------------------")
                    skip = switch_fragmentsize(gene, idx, OLS)
                    if skip:
                        continue
                    # Quality Control for overhangs from the same gene
                    # check_overhangs(gene, OLS)
                    SPINEgene.barcodeF.extend(compileF)
                    SPINEgene.barcodeR.extend(compileR)
                    compileF = []
                    compileR = []
                    idx = 0
                    continue
                # Store
                gene.genePrimer.append(SeqRecord(reverse, id=gene.geneid + "_geneP_Mut-" + str(idx + 1) + "_R",
                                                 description="Frag" + fragstart + "-" + fragend + ' ' + str(tmR) + 'C'))
                gene.genePrimer.append(SeqRecord(forward, id=gene.geneid + "_geneP_Mut-" + str(idx + 1) + "_F",
                                                 description="Frag" + fragstart + "-" + fragend + ' ' + str(tmF) + 'C'))
                # Determine missing double mutations
                beginning = int((frag[0] - SPINEgene.primerBuffer - sR) / 3)  # Region missing double mutations
                if beginning < 1:
                    beginning = 1
                end = ceil((frag[1] - SPINEgene.primerBuffer + sF) / 3)
                if end > ceil((gene.breaksites[-1] - SPINEgene.primerBuffer) / 3):
                    end = ceil((gene.breaksites[-1] - SPINEgene.primerBuffer) / 3)
                missingTmp = set()
                for site in range(beginning, end):  # Record these missing double mutations
                    for site2 in range(site + 1, end):
                        missingSites.append([site, site2])
                        missingTmp.add(site)
                        missingTmp.add(site2)
                        # missingTable[site][site2] = 0
                missingFragments.append([(frag[0] - 30) / 3, (frag[1] - 30) / 3, list(missingTmp)])
            if gene.unique_Frag[idx]:  # only for unique sequences
                # Create gene fragments with insertions
                tmF = 0
                tmR = 0
                count = 0
                tmpseq = gene.seq[frag[0] - 4-overlap : frag[1]+4+overlap].ungap('-') #4 is overhang for BsmBI
                offset = 4 + overlap
                tmpsequences = []
                # Create the mutations
                for i in range(offset, offset + frag[1] - frag[0], 3):
                    wt = [name for name, codon in gene.SynonymousCodons.items() if tmpseq[i:i + 3].upper() in codon]
                    for jk in (x for x in gene.aminoacids if x not in wt[0]):
                        p = [gene.usage[aa] for aa in gene.SynonymousCodons[jk]]  # Find probabilities
                        p = [xp if xp > 0.1 else 0 for xp in p]  # Remove probabilities below 0.1
                        p = [xp / sum(p) for xp in p]  # Normalize to 1
                        mutation = np.random.choice(gene.SynonymousCodons[jk], 1, p)  # Pick one codon
                        xfrag = tmpseq[0:i] + mutation[0] + tmpseq[i + 3:]  # Add mutation to fragment
                        # Check each cassette for more than 2 BsmBI and 2 BsaI sites
                        while xfrag.upper().count('GGTCTC') + xfrag.upper().count('GAGACC') > 2 | xfrag.upper().count('CGTCTC') + xfrag.upper().count('GAGACG') > 2:
                            print('Found BsaI and BsmBI sites')  # change codon
                            mutation = np.random.choice(gene.SynonymousCodons[jk], 1, p)  # Pick one codon
                            xfrag = tmpseq[0:i] + mutation + tmpseq[i + 3:]
                        tmpsequences.append(SeqRecord(xfrag,
                                                      id=gene.geneid + "_Mut" + fragstart + "-" + fragend + "_" + wt[0] + str(int((frag[0] + i + 3 - offset - SPINEgene.primerBuffer) / 3)) + jk,
                                                      description=''))

                if gene.num_frag_per_oligo > 1:
                    tmpsequences = combine_fragments(tmpsequences, gene.num_frag_per_oligo, gene.split)

                # add on barcodes
                tmpseq = tmpsequences[0].seq
                while tmF < SPINEgene.primerTm[0] or tmR < SPINEgene.primerTm[0]:  # swap out barcode if tm is low
                    difference = (SPINEgene.synth_len - (len(tmpseq) + 14))  # 14 bases is the length of the restriction sites with overhangs (7 bases each)
                    barF = SPINEgene.barcodeF.pop(0)
                    barR = SPINEgene.barcodeR.pop(0)
                    count += 1  # How many barcodes used
                    compileF.append(barF)
                    compileR.append(barR)
                    while difference / 2 > len(barF):
                        tmpF = SPINEgene.barcodeF.pop(0)
                        tmpR = SPINEgene.barcodeR.pop(0)
                        compileF.append(tmpF)
                        compileR.append(tmpR)
                        barF += tmpF
                        barR += tmpR
                        count += 1  # How many barcodes used
                    tmpfrag_1 = barF.seq[0:int(difference / 2)] + "CGTCTCC" + tmpseq[0:4]
                    tmpfrag_2 = tmpseq[-4:] + "GGAGACG" + barR.seq.reverse_complement()[0:difference - int(difference / 2)]
                    # primers for amplifying subpools
                    offset = int(difference / 2) + 11  # add 11 bases for type 2 restriction
                    primerF, tmF = find_fragment_primer(tmpfrag_1, offset)
                    primerR, tmR = find_fragment_primer(tmpfrag_2.reverse_complement(), (difference - offset + 22))
                group_oligos = []
                for sequence in tmpsequences:
                    combined_sequence = tmpfrag_1 + sequence.seq[4:-4] + tmpfrag_2
                    if gene.doublefrag == 0:
                        gene.oligos.append(SeqRecord(combined_sequence, id=sequence.id, description=''))
                    else:
                        grouped_oligos.append(SeqRecord(combined_sequence, id=sequence.id, description=''))

                # Store primers for gene fragment
                gene.barPrimer.append(SeqRecord(primerF, id=gene.geneid + "_oligoP_Mut-" + str(idx + 1) + "_F",
                                                description="Frag" + fragstart + "-" + fragend + "_" + str(tmF) + 'C'))
                gene.barPrimer.append(SeqRecord(primerR, id=gene.geneid + "_oligoP_Mut-" + str(idx + 1) + "_R",
                                                description="Frag" + fragstart + "-" + fragend + "_" + str(tmR) + 'C'))
                print('Barcodes used:' + str(count))
                print('Barcodes Remaining:' + str(len(SPINEgene.barcodeF)))
            if gene.doublefrag == 1:
                all_grouped_oligos.append(grouped_oligos)
            idx += 1
        # Resolve Double Fragment
        if gene.doublefrag == 1:
            while len(all_grouped_oligos) > 1:
                listOne = all_grouped_oligos.pop(0)
                listTwo = all_grouped_oligos.pop(0)
                while listOne and listTwo:
                    one = listOne.pop(0)
                    two = listTwo.pop(0)
                    combined_sequence = one.seq + two.seq.reverse_complement()
                    combined_id = one.id + two.id
                    gene.oligos.append(SeqRecord(combined_sequence, id=combined_id, description=''))
                if listOne or listTwo:
                    if listOne:
                        sequence = listOne.pop(0)
                    if listTwo:
                        sequence = listTwo.pop(0)
                    combined_id = sequence.id
                    combined_sequence = sequence.seq
                    difference = 230 - len(combined_sequence)
                    # print(len(tmpseq))
                    barF2 = SPINEgene.barcodeF.pop(0)
                    barR2 = SPINEgene.barcodeR.pop(0)
                    while difference / 2 > len(barF2):
                        barF2 += SPINEgene.barcodeF.pop(0)
                        barR2 += SPINEgene.barcodeR.pop(0)
                    combined_sequence2 = barF2.seq[0:int(difference / 2)] + combined_sequence + barR2.seq.reverse_complement()[0:difference - int(difference / 2)]
                    gene.oligos.append(SeqRecord(combined_sequence2, id=combined_id, description=''))
            if all_grouped_oligos:
                one = all_grouped_oligos
                while one:
                    sequence_one = one.pop(0)
                    combined_id = sequence_one.id
                    combined_sequence = sequence_one.seq
                    difference = 230 - len(combined_sequence)
                    # print(len(tmpseq))
                    barF2 = SPINEgene.barcodeF.pop(0)
                    barR2 = SPINEgene.barcodeR.pop(0)
                    while difference / 2 > len(barF2):
                        barF2 += SPINEgene.barcodeF.pop(0)
                        barR2 += SPINEgene.barcodeR.pop(0)
                    combined_sequence2 = barF2.seq[0:int(difference / 2)] + combined_sequence + barR2.seq.reverse_complement()[0:difference - int(difference / 2)]
                    gene.oligos.append(SeqRecord(combined_sequence2, id=combined_id, description=''))
        # Export files (fasta)
        # Missing Mutation Pairs
        import csv
        with open(os.path.join(folder.replace('\\', ''), gene.geneid + "_missing2Mutations.csv"), 'w') as csvfile:
            mutationwriter = csv.writer(csvfile, delimiter=',')
            mutationwriter.writerows(missingSites)
            mutationwriter.writerow('Fragment Info')
            mutationwriter.writerows(missingFragments)
        # Print table?
        # from tabulate import tabulate
        # print('Missing Double Mutation Table:')
        # print(tabulate(missingTable))
        # Fragments
        SeqIO.write(gene.oligos, os.path.join(folder.replace('\\', ''),
                                                 gene.geneid + "_Mut_Oligos.fasta"), "fasta")
        # Barcode Primers
        SeqIO.write(gene.barPrimer, os.path.join(folder.replace('\\', ''),
                                                 gene.geneid + "_Mut_Oligo_Primers.fasta"), "fasta")
        # Amplification Primers
        SeqIO.write(gene.genePrimer, os.path.join(folder.replace('\\', ''),
                                                  gene.geneid + "_Mut_Gene_Primers.fasta"), "fasta")
        # Record finished gene for aligned genes
        finishedGenes.extend([ii])


def combine_fragments(tandem, num_frag_per_oligo, split):
    tandem_seq = []
    barcodes = []
    if split:
        tmpF = SPINEgene.barcodeF.pop(0)
        tmpR = SPINEgene.barcodeR.pop(0)
    direction = -1
    while len(tandem) > num_frag_per_oligo:
        tmp = tandem.pop(0)
        tmp_tandem = tmp.seq
        tandem_id = tmp.id
        for x in range(num_frag_per_oligo - 1):
            if split:
                name = tmp.id
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += "GGAGACG" + tmpR + tmpF + "CGTCTCC" + tmp.seq  # concatenate and add cut sites with buffer
                    direction = -1
                else:
                    tmp_tandem += "GGAGACG" + tmpR + tmpF + "CGTCTCC" + tmp.seq.reverse_complement()  # concatenate and add cut sites with buffer
                    direction = 1
                tandem_id += '+' + tmp.id
                barcodes.append(SeqRecord(tmpR, id=name))
                barcodes.append(SeqRecord(tmpF, id=tmp.id))
            else:
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq  # concatenate and add cut sites with buffer
                    direction = -1
                else:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq.reverse_complement()
                    direction = 1
                tandem_id += '+' + tmp.id
        tandem_seq.append(SeqRecord(tmp_tandem, id=tandem_id, description=''))

    if tandem:
        direction = -1
        print(len(tandem))
        tmp = tandem.pop(0)
        tmp_tandem = tmp.seq
        tandem_id = tmp.id
        while tandem:
            if split:
                name = tmp.id
                tmp = tandem.pop(0)
                if direction == 1:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq  # concatenate and add cut sites with buffer
                    direction = 1
                else:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq.reverse_complement()  # concatenate and add cut sites with buffer
                    direction = -1
                tandem_id += '+' + tmp.id
                barcodes.append(SeqRecord(tmpR, id=name))
                barcodes.append(SeqRecord(tmpF, id=tmp.id))
            else:
                tmp = tandem.pop(0)
                if direction == -1:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq  # concatenate and add cut sites with buffer
                    direction = 1
                else:
                    tmp_tandem += "GGAGACG" + "ACGT" + "CGTCTCC" + tmp.seq
                    direction = -1
                tandem_id += '+' + tmp.id
        difference = len(tandem_seq[-1].seq) - len(tmp_tandem)
        barF = SPINEgene.barcodeF.pop(0)
        barR = SPINEgene.barcodeR.pop(0)
        while difference / 2 > len(barF):
            barF += SPINEgene.barcodeF.pop(0)
            barR += SPINEgene.barcodeR.pop(0)
        tmpfrag = barF.seq[0:int(difference / 2)] + tmp_tandem + barR.seq.reverse_complement()[0:difference - int(difference / 2)]
        tandem_seq.append(SeqRecord(tmpfrag, id=tandem_id, description=''))
        print('Partial sequence' + str(len(tmpfrag)))
    return tandem_seq


def print_all(OLS, folder=''):
    if not isinstance(OLS[0], SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    alloligos = []
    allprimers = []
    for obj in OLS:
        try:
            alloligos.extend(obj.oligos)
            allprimers.extend(obj.barPrimer)
            allprimers.extend(obj.genePrimer)
        except AttributeError:
            print(obj.geneid + " has not been processed")
    # Remove redundant sequences?
    SeqIO.write(alloligos, os.path.join(folder.replace('\\', ''), "All_Oligos.fasta"), "fasta")
    SeqIO.write(allprimers, os.path.join(folder.replace('\\', ''), "All_Primers.fasta"), "fasta")


def post_qc(OLS):
    if not isinstance(OLS[0], SPINEgene):
        raise TypeError('Not an instance of the SPINEgene class')
    # Post QC
    # all_oligos = list(SeqIO.parse(r'/Users/Luna/Google Drive/BIOLOGYPLUS/DNA DATABASE/P03 - SWITCH INSERTION/OLS_Library_Generation/OLS_003/All_Oligos.fasta','fasta'))
    all_oligos = []
    # all_barPrimers = list(SeqIO.parse(r'/Users/Luna/Google Drive/BIOLOGYPLUS/DNA DATABASE/P03 - SWITCH INSERTION/OLS_Library_Generation/OLS_003/All_Primers.fasta','fasta'))
    all_barPrimers = []
    for obj in OLS:
        try:
            all_oligos.extend(obj.oligos)
            all_barPrimers.extend(obj.barPrimer)
        except AttributeError:
            print(obj.geneid + " has not been processed")

    print("Running QC for barcode primer specificiy")
    cassetteSet = set(all_oligos[0].id[:-6])
    uCassette = [SeqRecord(all_oligos[0].seq, id=all_oligos[0].id[:-6])]
    for idx in range(len(all_oligos)):
        if not all_oligos[idx].id[:-6] in cassetteSet:
            uCassette.append(SeqRecord(all_oligos[idx].seq, id=all_oligos[idx].id[:-6]))
        cassetteSet.add(all_oligos[idx].id[:-6])
    grouped = iter(all_barPrimers)
    grouped = zip(grouped, grouped)  # create the combinatorial comparisons
    nonspecific = []
    # iterate over every barcode primer pair and match to each oligo to check for nonspecific amplification
    for idxPrime, primers in enumerate(grouped):  # iterate over every barcode primer pair
        print("Checking primer set:" + primers[0].id[:-2])
        for idxCassette, fragment in enumerate(uCassette):  # iterate over every OLS oligo
            if primers[0].id.split('_')[2] != all_oligos[idxCassette].id.split('_')[2]:  # ignore designed annealing (same name)
                fragname = fragment.id
                fragment = fragment.seq
                non = [[False], [False]]
                for idxDirection, primer in enumerate(primers):
                    primername = primer.id
                    primer = primer.seq
                    for i in range(
                            len(fragment) - len(primer)):  # iterate over the length of the oligo for a binding site
                        match = [primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))]
                        first = 10
                        for k in range(len(match) - 3):
                            if (match[k] and match[k + 1] and match[k + 3]) or (
                                    match[k] and match[k + 1] and match[k + 2]):
                                first = k
                                break
                        if sum(match[first:]) > len(primer[first:]) * 0.8 and sum(match[first:]) > 6 and match[-1]:  # string compare - sum of matched nt is greater than 80%
                            try:
                                melt = mt.Tm_NN(primer[first:], c_seq=fragment[i + first:i + len(primer)].complement(),
                                                nn_table=mt.DNA_NN2, de_table=mt.DNA_DE1, imm_table=mt.DNA_IMM1)
                                # if melt>20:
                                #     print('Found non-specfic match:'+fragment[i:i+len(primer)])
                                #     print('                 primer:'+primer+' Tm:'+str(round(melt,1)))
                                if melt > 35:
                                    non[0].append(True)
                            except ValueError as valerr:
                                # print(str(valerr)+". Please check position manually:"+str(i+1)+" forward")
                                # print('Primer:'+primer)
                                # print('Match: '+fragment[i:i+len(primer)])
                                pass
                    fragment = fragment.reverse_complement()
                    for i in range(len(fragment) - len(primer)):
                        match = [primer[j].lower() == fragment[i + j].lower() for j in range(len(primer))]
                        first = 10
                        for k in range(0, len(match) - 3, 1):
                            if match[k] and match[k + 1] and match[k + 3]:
                                first = k
                                break
                        if sum(match[first:]) > len(primer[first:]) * 0.8 and sum(match[first:]) > 6 and match[-1]:  # string compare - sum of matched nt is greater than 80%
                            try:
                                melt = mt.Tm_NN(primer[first:], c_seq=fragment[i + first:i + len(primer)].complement(),
                                                nn_table=mt.DNA_NN2, de_table=mt.DNA_DE1, imm_table=mt.DNA_IMM1)
                                # if melt > 20:
                                #     print('Found non-specfic match:'+fragment[i:i+len(primer)])
                                #     print('                 primer:'+primer+' Tm:'+str(melt))
                                if melt > 35:
                                    non[1].append(True)
                            except ValueError as valerr:
                                # print(str(valerr)+". Please check position manually:"+str(i+1)+" reverse")
                                # print('Primer:'+primer)
                                # print('Match: '+fragment[i:i+len(primer)])
                                pass
                    if sum(non[0]) == 0 and sum(non[1]) == 0:
                        break
                    if sum(non[0]) > 0 and sum(non[1]) > 0:
                        nonspecific.append([primername, fragname])
                        print("Found Non-specific Amplification")
    if nonspecific:
        print('Nonspecific Primers: (Manually changing primer sequence recommended)')
        print(nonspecific)
    else:
        print('No non-specific primers detected')
