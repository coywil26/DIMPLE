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