# Q1 II: Reading from a fasta file using the SeqIO
# module in biopython

# from Bio import SeqIO
#
# file_name = 'ATdreb2a.fasta'
# seq_object = SeqIO.read(file_name, "fasta")
#
# print(f"Sequence ID: {seq_object.id}\n"
#       f"Description: {seq_object.description}\n"
#       f"Sequence:\n"
#       f"{seq_object.seq}\n"
#       f"Sequence length: {len(seq_object.seq)}")

# Q1 III: Carrying out a web blast search using
# NCBIWWW qblast() function
# from Bio import SeqIO
# from Bio.Blast import NCBIWWW
#
# file_name = 'ATdreb2a.fasta'
# seq_object = SeqIO.read(file_name, "fasta")
#
# print("Retrieving blast result...")
# result = NCBIWWW.qblast("blastn", "nt", seq_object.format("fasta"))
# print("Blast result retrieved successfully.")
#
# with open("dreb2a_blast.xml", "w") as output:
#       output.write(result.read())
#
# result.close()
# print("Blast result stored in output file successfully.")


# Q1 IV: Extracting the details related to the hit sequences
# from a xml file containing the blast search results using NCBIXML
# module in biopython
# from Bio.Blast import NCBIXML
#
# results = open("dreb2a_blast.xml")
#
# blast_record = NCBIXML.read(results)
#
# evalue_thresh = 0.05
#
# for hit in blast_record.alignments:
#
#       print(f"Blast hit title: {hit.title}\n"
#             f"Hit sequence length: {hit.length}\n")
#
#       hsp_num = 0
#       for hsp in hit.hsps:
#             hsp_num += 1
#             if hsp.expect < evalue_thresh:
#                   print(f"HSP number: {hsp_num}\n"
#                         f"Alignment length: {len(hsp.sbjct)}\n"
#                         f"e value: {hsp.expect} \n"
#                         f"Score: {hsp.score}\n"
#                         f"Subject sequence: \n{hsp.sbjct}\n")

# Q1 V: Using regular expressions to identify ABRE elements within the
# alignment/ hit sequences in the blast search results
import re
from Bio.Blast import NCBIXML

# Method to find ABRE elements in a given sequence
# def find_ABRE(hsp, seq):
#       pattern = re.compile(r"((?:C|T)ACGT(?:G|T)C)")
#
#       mo_1 = re.finditer(pattern, seq)
#
#       # Dictionary that stores the relevant matching ABRE element
#       # and the spanning length of that ABRE element in the sequence
#       # in each if the hsp of a single hit sequence
#       ABRE_ele = {}
#
#       # Finding the exact location of the elements in the hit sequence
#       # and appending it to ABRE_ele dictionary
#       for match in mo_1:
#             if hsp.frame == 1:
#                   ABRE_ele[match.group()] = (match.span()[0] + hsp.sbjct_start, match.span()[1] + hsp.sbjct_start)
#             if hsp.frame == -1:
#                   ABRE_ele[match.group()] = (match.span()[0] + hsp.sbjct_end, match.span()[1] + hsp.sbjct_end)
#       return ABRE_ele


# Finding the ABRE elements in the hit sequences obtained from our previous
# blast search
# results = open("dreb2a_blast.xml")
# blast_record = NCBIXML.read(results)
#
# hits_with_ABRE = {}
# hits = 0
# for hit in blast_record.alignments:
#
#       ABRE_in_hit = {}
#
#       hits += 1
#
#       hsp_num = 0
#       for hsp in hit.hsps:
#
#             hsp_num += 1
#             if find_ABRE(hsp, hsp.sbjct) != {}:
#                   ABRE_in_hit[f'HSP {hsp_num}'] = find_ABRE(hsp, hsp.sbjct)
#
#       acc = hit.title.split('|')[3]
#       if ABRE_in_hit != {}:
#             hits_with_ABRE[acc] = ABRE_in_hit
#
# print(f"\nHit sequences wit ABRE element: \n")
#
# for acc, value in hits_with_ABRE.items():
#       print(f"{acc}: {value}")
#
# print(f"\nTotal number of hits: {hits}\n"
#       f"Number of hits with ABRE element: {len(hits_with_ABRE)}")










