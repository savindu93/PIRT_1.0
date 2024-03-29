import re

pattern1 = r"[^AGTC\n]"
pattern2 = "(XM)|(NM)"

f = open("OSDREB1A.txt", 'r')

mRNA = ""

# Identifying the mRNA sequences from the amino acid sequence using regular
# expression and the patterns defined above
for line in f:
    match = re.search(pattern1, line)
    if match2 := re.search(pattern2, line):
        mRNA += line.strip() + "| Transcribed\n"
    if not match:
        line = line.replace("T", "U")
        mRNA += line.strip()

print(mRNA)

f1 = open("OSDREB1A_mRNA_1.fasta", 'w')
f1.write(mRNA)
f.close()
f1.close()
