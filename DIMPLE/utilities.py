from Bio import pairwise2

def find_mutations(mutations, wt_seq, typeII_RE, typeII_RE_R, RE_gap):
    mutation_file = r"C:\Users\test.fasta"
    typeII_RE = 'CGTCTC'  # I am using BsmbI
    typeII_RE_R = 'GAGACG'
    wt_seq = 'ATG'  # coding sequencing only
    mutations = []
    RE_gap = 6
    with open(mutation_file, 'r') as file:
        row = next(file)
        name = row.split('_')[6].strip()
        mutations.append('>'+name)
        for row in file:
            full_row = ''
            while '>' != row[0]:
                full_row += row.strip()
                try:
                    row = next(file)
                except:
                    break
            if full_row == '':
                break
            #codon = ''.join([x for x in name if x.isdigit()])
            tmp_split = full_row.strip().split(typeII_RE)[1][RE_gap:]  # minimum of 5 bases
            row_split = tmp_split.split(typeII_RE_R)[0][:-RE_gap]
            s = pairwise2.align.localms(row_split, wt_seq, 2, 0, -10, -10)[0][0]
            start = len(s)-len(s.lstrip('-'))
            remainder = start % 3
            oligo_start = remainder
            for codon in range(oligo_start, len(row_split), 3):
                if row_split[codon:codon+3] != wt_seq[start+codon:start+codon+3]:
                    mutations.append(row_split[codon:codon+3])
            if '>' in row:
                name = row.split('_')[6].strip()
                mutations.append('>'+name)

    with open(r"C:\Users\test_mutations.csv", 'w') as file:
        for mut in mutations:
            file.write(mut+'\n')

    return mutations


def parse_custom_mutations(mutation_text):
    custom_mutations = {}
    for set in mutation_text:
        set = set.split(':')
        if set[1] == 'All':
            if '-' in set[0]:
                for i in range(int(set[0].split('-')[0]), int(set[0].split('-')[1]) + 1):
                    custom_mutations[i] = 'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y'
            else:
                custom_mutations[int(set[0])] = 'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y'
        else:
            if '-' in set[0]:
                for i in range(int(set[0].split('-')[0]), int(set[0].split('-')[1]) + 1):
                    custom_mutations[i] = set[1]
            else:
                # if mutation exists, add to it
                if int(set[0]) in custom_mutations.keys():
                    custom_mutations[int(set[0])] = custom_mutations[int(set[0])] + ',' + set[1]
                else:
                    custom_mutations[int(set[0])] = set[1]
    return custom_mutations


def codon_usage(usage):
    if usage == 'ecoli':
        usage_table = {
            'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.09,
            'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
            'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
            'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
            'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGA': 0.3, 'TGG': 1,
            'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
            'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
            'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15
        }  # E.coli codon usage table
    elif usage == 'human':
        usage_table = {
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
        usage_table = usage
    return usage_table

