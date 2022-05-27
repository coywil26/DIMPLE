from Bio import pairwise2
mutation_file = r"C:\Users\david.nedrud\Documents\GitHub\SPINE_DN\tests\DHFR_DMS_Oligos.fasta"
typeII_RE = 'CGTCTC' # I am using BsmbI
typeII_RE_R = 'GAGACG'
wt_seq = 'ATGGTGCGTCCGTTAAACTGCATTGTCGCTGTGTCACAAAACATGGGGATTGGGAAGAATGGCGACCTCCCATGGCCGCCGCTGCGGAACGAATTCAAATACTTCCAACGGATGACGACGACCTCGTCAGTGGAGGGCAAGCAAAACCTCGTGATTATGGGCCGCAAGACGTGGTTCAGTATCCCGGAAAAGAACCGCCCTCTTAAAGACCGGATCAATATTGTGCTTTCCCGGGAGTTGAAAGAACCGCCGCGTGGCGCACATTTCCTGGCAAAAAGTCTTGACGACGCACTTCGTTTAATTGAGCAACCGGAGTTGGCATCAAAAGTGGACATGGTGTGGATTGTTGGCGGCAGTAGTGTTTACCAGGAAGCCATGAACCAGCCAGGCCATCTTCGCTTGTTTGTTACACGGATTATGCAAGAGTTTGAGAGCGACACATTTTTCCCGGAGATCGACCTGGGCAAGTACAAACTTCTCCCTGAATATCCTGGGGTGTTATCAGAAGTCCAAGAAGAGAAAGGGATCAAGTACAAGTTTGAAGTTTACGAAAAGAAGGATTAA'  # coding sequencing only
mutations = []

with open(mutation_file,'r') as file:
    row = next(file)
    name = row.split('_')[2].strip()
    mutations.append('>'+name)
    for row in file:
        full_row = ''
        while '>' != row[0]:
            full_row += row.strip()
            row = next(file)
        if full_row == '':
            break
        #codon = ''.join([x for x in name if x.isdigit()])
        tmp_split = full_row.strip().split(typeII_RE)[1][5:]  # minimum of 5 bases
        row_split = tmp_split.split(typeII_RE_R)[0][:-5]
        s = pairwise2.align.localms(row_split, wt_seq,2,0,-10,-10)[0][0]
        start = len(s)-len(s.lstrip('-'))
        remainder = start%3
        oligo_start = remainder
        for codon in range(oligo_start,len(row_split),3):
            if row_split[codon:codon+3] != wt_seq[start+codon:start+codon+3]:
                mutations.append(row_split[codon:codon+3])
        name = row.split('_')[2].strip()
        mutations.append('>'+name)

with open("geneA_mutations.csv",'w') as file:
    for mut in mutations:
        file.write(mut+'\n')
