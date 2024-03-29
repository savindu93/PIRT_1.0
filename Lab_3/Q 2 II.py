#Q 2 II
f3 = open("../Lab_4/codon_table.txt", 'r')

codon_dict = {}

# Storing the codons translating to its specific amino acid with its letter in a dicttionary
for line in f3:
    if line.startswith('#'):
        continue
    else:
        amino_acid_info = line.split()
        if len(codon_dict) == 0:
            codon_dict[amino_acid_info[2]] = [amino_acid_info[0]]
        else:
            if amino_acid_info[2] in codon_dict:
                codon_dict[amino_acid_info[2]].append(amino_acid_info[0])
            else:
                codon_dict[amino_acid_info[2]] = [amino_acid_info[0]]

print(f"Codon dictionary: \n", codon_dict)
print("\n")

f4 = open("OSDREB1A_mRNA.fasta", 'r')

mRNA_seq = ""
codons = []
amino_acid_seq = ""

# Create codons from the 1st base of the mRNA sequence
for line in f4:
    if not line.startswith('>'):
        mRNA_seq += line

i = 0
while i != len(mRNA_seq):
    codons.append(mRNA_seq[i: i+3])
    i += 3

print("Codons: \n", codons)
print("\n")

#Translating the mRNA seq to an amino acid sequence
j = 0
while (j != len(codons)) and (codons[j] not in codon_dict['O']):
    amino_acid = [key for key, value in codon_dict.items() if codons[j] in value]
    amino_acid_seq += amino_acid[0]
    j += 1

print(f"Amino acid sequence: {amino_acid_seq}")
print("\n")
print(f"Length of amino acid sequence: {len(amino_acid_seq)}")
print("\n")

amino_acid_fasta = f"> Translated amino acid sequence\n {amino_acid_seq}"

f5 = open("amino_acid.fasta", 'w')
f5.write(amino_acid_fasta)
