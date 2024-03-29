# Q 1) IV
# f = open("OSDREB1A.txt", 'r')
#
# mRNA = ""
#
# for line in f:
#     if "XM_015755426.2" in line:
#             mRNA += line.strip() + "| Transcribed\n"
#     if len(mRNA) != 0 and not line.startswith('>') :
#         line = line.replace("T", "U")
#         mRNA += line.strip()
#
# print(mRNA)
#
# f1 = open("OSDREB1A_mRNA.fasta", 'w')
# f1.write(mRNA)

# f.close()
# f1.close()

### With regular expression

# import re
# pattern1 = r"[^AGTC\n]"
# pattern2 = "(XM)|(NM)"
#
# f = open("OSDREB1A.txt", 'r')
#
# non_mRNA = ""
# mRNA = ""
#
# for line in f:
#     match = re.search(pattern1, line)
#     if match2 := re.search(pattern2, line):
#         mRNA += line.strip() + "| Transcribed\n"
#     if not match:
#         line = line.replace("T", "U")
#         mRNA += line.strip()
#
# print(mRNA)
# print(non_mRNA)
#
# f1 = open("OSDREB1A_mRNA_1.fasta", 'w')
# f1.write(mRNA)
#
# f.close()
# f1.close()

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

f4 = open("OSDREB1A_mRNA.fasta", 'r')

mRNA_seq = ""
codons = []
amino_acid_seq = ""

# Create codons from 1st base
for line in f4:
    if not line.startswith('>'):
        mRNA_seq += line

i = 0
while i != len(mRNA_seq):
    codons.append(mRNA_seq[i: i+3])
    i += 3

print("Codons: \n", codons)

# creating all possible reading frames
# j = 0
# k = 0
# while j < 3:
#     i = j
#     print(f"Reading frame {j+1}")
#     while i < len(mRNA_seq):
#         if mRNA_seq[i: i+3] == 'AUG':
#             k = i
#             reading_frame = []
#             while mRNA_seq[k: k+3] not in ['UAA', 'UAG', 'UGA']:
#                 if k == len(mRNA_seq) - 3:
#                     reading_frame = []
#                     break
#                 reading_frame.append(mRNA_seq[k: k+3])
#                 k += 3
#             print(f"seq: {reading_frame}")
#             if len(reading_frame) != 0:
#                 codons.append(reading_frame)
#         i += 3
#     j += 1

#print(codons)

#Translating the mRNA seq to an amino acid sequence
j = 0
while (j != len(codons)) and (codons[j] not in codon_dict['O']):
    amino_acid = [key for key, value in codon_dict.items() if codons[j] in value]
    amino_acid_seq += amino_acid[0]
    j += 1

print(f"Amino acid sequence: {amino_acid_seq}")

f5 = open("amino_acid.fasta", 'w')
f5.write(amino_acid_seq)



