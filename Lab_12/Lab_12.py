# Q1) I
# accs = ["AAK43967.1", "AED90870.1", "NP_567720.1","AAK59861.1"]
#
# for acc in accs:
#
#     stream = Entrez.efetch(db = 'protein', id = acc, rettype = 'gb', retmode = 'text')
#     record = SeqIO.read(stream, 'genbank')
#     stream.close()
#
#     seq = ""
#     seq += f">{record.id}\n" \
#            f"{record.seq}"
#
#     file = open(f"{acc}.fasta", 'w')
#     file.write(seq)
#     file.close()

# import re
# from Bio import SeqIO
#
# pattern = re.compile("(WGKWV)|(AAEIR)")
#
# files = ["AAK43967.1", "AED90870.1", "NP_567720.1","AAK59861.1"]
#
#
# for file in files:
#
#     with open(f"{file}.fasta",'r') as handle:
#
#         seq_o = SeqIO.read(handle,"fasta")
#         seq = str(seq_o.seq)
#         mo = re.search(pattern, seq)
#
#         if mo:
#             out_file = open("AP2_advanced_headers.txt",'a')
#             out_file.write(f"{seq_o.id}\n")
#             out_file.close()


# Cds_seq_retrieve
from Bio import Entrez
from Bio import SeqIO
import sys

acc = "NM_000188.3"

stream = Entrez.efetch(db = 'nucleotide', id = acc, rettype = 'gb', retmode = 'text')
record = SeqIO.read(stream, 'genbank')
stream.close()

seq = ""
seq += f">{record.id}\n" \
       f"{record.seq}"

print(seq)

cds_feature = [feature for feature in record.features if feature.type == "CDS"]

if cds_feature:

       # Extract the coding sequence from the CDS feature
       cds_sequence = cds_feature[0].extract(record.seq)
       print(cds_sequence)


# file = open(f"{acc}_cds_seq.fasta", 'w')
# file.write(seq)
# file.close()

# # Transcribe sequence
# with open(file,'r') as file:
#
#        dna = SeqIO.read(file).seq
#        transcribed_dna = dna.replace("T","U")
#
#        seq = f">{dna.id}-transcribed\n" \
#              f"{transcribed_dna}\n"
#
#        mRNA_file = open("mRNA_seq.fasta",'w')
#        mRNA_file.write(seq)
#        mRNA_file.close()
#
#
# # Translate sequence
# with open(file,'r') as file:
#
#        mRNA = SeqIO.read(file).seq
#        translated = mRNA.translate()
#
#        seq = f">{mRNA.id}-translated\n" \
#              f"{translated}\n"
#
#        aa_file = open("aa_seq.fasta",'w')
#        aa_file.write(seq)
#        aa_file.close()
#
# # Analyze aa sequence
# import sys
# from Bio.SeqUtils import molecular_weight
#
# file = sys.argv[1:]
#
# with open(file,'r') as file:
#
#        seq = SeqIO.read(file).seq
#
#        # Calculate Total molecular weight (mw) of the amino acid (aa) sequence/
#        # Calculate the alanine and glycine percentage
#
#        tot_length = len(seq) # total length of aa sequence
#
#        tot_mw =  molecular_weight(seq) # total mw
#
#        A = 0  # total alanine count
#        G = 0  # total glycine count
#        for aa in seq:
#            if aa is 'A': A+= 1
#            if aa is 'G': G += 1
#
#        A_content = (A/tot_length) * 100 # total alanine percentage
#        G_content = (G/tot_length) * 100 # total glycine percentage
#
#        stats = f"Amino acid;\n" \
#                f"Length: {tot_length}\n" \
#                f"Molecular weight: {tot_mw}\n" \
#                f"Alanine percentage: {A_content}\n" \
#                f"Glycine percentage: {G_content}\n"
#
#        file = open("aa_stats.txt",'w')
#        file.write(stats)
#        file.close()


























