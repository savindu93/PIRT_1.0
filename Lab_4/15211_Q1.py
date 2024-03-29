# Q1 III; Calculate AT content in a sequence
def cal_AT(dna_seq):

        A = 0
        T = 0

        for base in dna_seq:
            if 'A' in base:
                A += 1
            if 'T' in base:
                T += 1

        return f"AT content: {(A + T)/len(dna_seq)}"

# Q1 IV; Split multiple fasta sequences and return a dictionary
def extract_seq(fasta_file):

        seq_dict = {}
        key = ""
        value = ""

        # Identifying the end of file in terms of
        # number of lines
        eof = 0
        with open(fasta_file, 'r') as file:

            for line in file:
                eof += 1

        with open(fasta_file, 'r') as file:

            line_no = 0
            for line in file:
                line_no += 1
                if line.strip() == "" or line_no == eof:
                    if line_no == eof:
                        value += line.strip()

                    # Adding a new key value pair in the sequence
                    # dictionary each time an empty line or the
                    # end of file is reached in the fasta file

                    seq_dict[key] = value
                    key = ""
                    value = ""
                elif line.startswith('>'):
                    key += line.strip('>').strip()
                else:
                    value += line.strip()

        return seq_dict

# Q1 V; Check whether a sequence is a protein, DNA or RNA
def check_seq_type(seq):
        import re

        pattern_1 = "[^AGTC]" # checks for protein or DNA sequence
        pattern_2 = "[^AGUC]" # checks for a mRNA sequence

        if re.search(pattern_1, seq) and re.search(pattern_2, seq):
            return "Protein"
        elif not re.search(pattern_1, seq):
            return "DNA"
        elif not re.search(pattern_2, seq):
            return "mRNA"

def print_dict(dict):

        for k, v in dict.items():
            print(f"{k} : {v}")

def print_dict_1(dict):

       for k, v in dict.items():
           print(k)
           print(f"Sequence type: {check_seq_type(v)}")
           if "DNA" in check_seq_type(v):
               print(cal_AT(v))


# Q1 III
dna_seq = "AATGCAATCGC"
print(cal_AT(dna_seq))

# Q1 VI
seq_dict = extract_seq('OSDREB_sequences_1.fasta')
print_dict_1(seq_dict)

