class Sequence:


    seq_count = 0

    # Constructor method to creates an instance of sequence class
    # or a Sequence object
    def __init__(self, gene_ID, gene_name, seq_type, seq_length, sp_name,
                 sbsp_name):

        self.gene_ID = gene_ID
        self.gene_name = gene_name
        self.seq_type = seq_type
        self.seq_length = seq_length
        self.sp_name = sp_name
        self.sbsp_name = sbsp_name
        Sequence.seq_count += 1

    # Q2 ix
    @staticmethod
    def fasta_split(fasta_file):

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

                    list = key.split('-')

                    # Adding a new key value pair in the sequence
                    # dictionary each time an empty line or the
                    # end of file is reached in the fasta file

                    seq_dict[list[0]] = {'gene_name': list[0],
                                         'gene_ID': list[1],
                                         'sp_name': list[2],
                                         'sbsp_name': list[3],
                                         'seq': value}

                    # Emptying the key, value variables to store
                    # next key and value

                    key = ""
                    value = ""
                elif line.startswith('>'):
                    key += line.strip('>').strip()
                else:
                    value += line.strip()

        return seq_dict

    # Q2 xi
    def get_character_count(self, seq):

        character_dict = {}

        for base in seq:
            if 'A' in base:
                if 'A' not in character_dict:
                    character_dict.update({'A': 1})
                else:
                    character_dict['A'] += 1
            if 'G' in base:
                if 'G' not in character_dict:
                    character_dict.update({'G': 1})
                else:
                    character_dict['G'] += 1
            if 'C' in base:
                if 'C' not in character_dict:
                    character_dict.update({'C': 1})
                else:
                    character_dict['C'] += 1
            if 'T' in base:
                if 'T' not in character_dict:
                    character_dict.update({'T': 1})
                else:
                    character_dict['T'] += 1

        return character_dict


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
                    seq_dict[key] = value
                    key = ""
                    value = ""
                elif line.startswith('>'):
                    key += line.strip('>').strip()
                else:
                    value += line.strip()

        return seq_dict

    @staticmethod
    # Q1 V; Check whether a sequence is a protein, DNA or RNA
    def check_seq_type(seq):
        import re

        pattern_1 = "[^AGTC]" # check sequence is a protein or DNA
        pattern_2 = "[^AGUC]" # check sequence is a mRNA

        if re.search(pattern_1, seq) and re.search(pattern_2, seq):
            return "Protein"
        elif not re.search(pattern_1, seq):
            return "DNA"
        elif not re.search(pattern_2, seq):
            return "mRNA"

    @classmethod
    # Function to print a dictionary containing
    # sequence data
    def print_seq_dict(cls, dict):

        for k, v in dict.items():
            print(k)
            print(cls.check_seq_type(v['seq']))
            if "DNA" in cls.check_seq_type(v['seq']):
                print(cls.cal_AT(v['seq']))

    # Creates a sequence object by taking a name of sequence
    # in this case gene and the sequence dictionary that contains
    # the information related to the sequence as args and calling
    # the Sequence() constructor method
    @staticmethod
    def create_seq_object(seq_dict, gene_name):

        return Sequence(seq_dict[gene_name]['gene_ID'],
                        seq_dict[gene_name]['gene_name'],
                        Sequence.check_seq_type(seq_dict[gene_name]['seq']),
                        len(seq_dict[gene_name]['seq']),
                        seq_dict[gene_name]['sp_name'],
                        seq_dict[gene_name]['sbsp_name'])


if __name__ == '__main__':

    seq_dict = Sequence.fasta_split('OSDREB_sequences_2.fasta')

    DREB1A = Sequence.create_seq_object(seq_dict, 'DREB1A_CDS')

    print(f"Gene_ID: {DREB1A.gene_ID} \n"
          f"Sequence length: {DREB1A.seq_length} \n"
          f"Sequence type: {DREB1A.seq_type}")

    gene_seq = seq_dict['DREB1A_CDS']['seq']
    print(f"\nBase count for DREB1A coding sequence: "
          f"\n{DREB1A.get_character_count(gene_seq)}")



