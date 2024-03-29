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

    # Method to split both protein and DNA seqeunces from fasta
    # file and storing it in one sequence dictionary
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

                    if 'CDS' in list[0]:

                        seq_dict[list[0]] = {'gene_name': list[0],
                                             'gene_ID': list[1],
                                             'sp_name': list[2],
                                             'sbsp_name': list[3],
                                             'seq': value}

                    if 'P' in list[0]:

                        seq_dict[list[0]] = {'protein_name': list[0],
                                             'gene_ID': list[1],
                                             'sp_name': list[2],
                                             'sbsp_name': list[3],
                                             'UniProt_ID':list[4],
                                             'reviewed_status':list[5],
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
    def get_ATcontent(seq):

        A = 0
        T = 0

        for base in seq:
            if 'A' in base:
                A += 1
            if 'T' in base:
                T += 1

        return f"{(A + T) / len(seq)}"

    @staticmethod
    # Q1 V; Check whether a sequence is a protein, DNA or RNA
    def check_seq_type(seq):
        import re

        pattern_1 = "[^AGTC]"  # check sequence is a protein or DNA
        pattern_2 = "[^AGUC]"  # check sequence is a mRNA

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


class DNAseq(Sequence):

    def __init__(self, gene_ID, gene_name, seq_type, seq_length, sp_name,
                 sbsp_name, AT_content, TraSc_seq):

        super().__init__(gene_ID, gene_name, seq_type, seq_length, sp_name,
                         sbsp_name)

        self.AT_content = AT_content
        self.TraSc_seq = TraSc_seq

    # The following method automates the creation of a DNA sequence
    # object rather than having to manually enter the sequence each
    # time. It uses the sequence information stored in seq_dict
    # dictionary and the relevant gene name and finds the information
    # related to that gene in the dictionary to create the sequence
    # object.
    @staticmethod
    def create_DNAseq_object(seq_dict, gene_name):

        DNA_seq = seq_dict[f"{gene_name}_CDS"]['seq']

        AT_content = DNAseq.get_ATcontent(DNA_seq)
        TraSc_seq =  DNAseq.transcribe_sequence(DNA_seq)

        return DNAseq(seq_dict[f"{gene_name}_CDS"]['gene_ID'],
                        seq_dict[f"{gene_name}_CDS"]['gene_name'],
                        Sequence.check_seq_type(DNA_seq),
                        len(DNA_seq),
                        seq_dict[f"{gene_name}_CDS"]['sp_name'],
                        seq_dict[f"{gene_name}_CDS"]['sbsp_name'],
                        AT_content,
                        TraSc_seq)

    def get_ATcontent(dna_seq):

        A = 0
        T = 0

        for base in dna_seq:
            if 'A' in base:
                A += 1
            if 'T' in base:
                T += 1

        return f"{(A + T) / len(dna_seq)}"

    def transcribe_sequence(seq):

        return seq.replace('T','U')

class mRNAseq(Sequence):

    __amino_acid_codons = {}

    def __init__(self, gene_ID, gene_name, seq_type, seq_length, sp_name,
                 sbsp_name, AT_content, TraLa_seq):

        super().__init__(gene_ID, gene_name, seq_type, seq_length, sp_name,
                         sbsp_name)

        self.AT_content = AT_content
        self.TraLa_seq = TraLa_seq

    # Method to get the AT content of a mRNA sequence
    def get_ATcontent(mRNA_seq):

        DNA_seq = mRNA_seq.replace('U','T')

        A = 0
        T = 0

        for base in DNA_seq:
            if 'A' in base:
                A += 1
            if 'T' in base:
                T += 1

        return f"{(A + T) / len(mRNA_seq)}"

    # Method to store codon infomration from an external file
    # into a dictionary
    @classmethod
    def upload_Codons(cls, file_name):

        codons_dict = cls.__amino_acid_codons

        f3 = open(file_name, 'r')

        # Storing the codons translating to its specific amino acid with its
        # letter in a dictionary
        for line in f3:
            if line.startswith('#'):
                continue
            else:
                amino_acid_info = line.split()
                if len(codons_dict) == 0:
                    codons_dict[amino_acid_info[2]] = [amino_acid_info[0]]
                else:
                    if amino_acid_info[2] in codons_dict:
                        codons_dict[amino_acid_info[2]].append(amino_acid_info[0])
                    else:
                        codons_dict[amino_acid_info[2]] = [amino_acid_info[0]]

        return codons_dict

    # The following method automates the creation of a mRNA sequence
    # object rather than having to manually enter the sequence each
    # time. It uses the sequence information stored in seq_dict
    # dictionary and the relevant gene name and finds the information
    # related to that gene in the dictionary to create the sequence
    # object.
    @staticmethod
    def create_mRNAseq_object(seq_dict, gene_name, seq = 'NA'):

        if seq != 'NA':
            if Sequence.check_seq_type(seq) == 'DNA':
                mRNA_seq = seq.replace('T', 'U')
            else:
                mRNA_seq = seq

        if seq == 'NA':
            mRNA_seq = seq_dict[f"{gene_name}_CDS"]['seq'].replace('T', 'U')

        AT_content = mRNAseq.get_ATcontent(mRNA_seq)
        TraLa_seq = mRNAseq.translate_sequence(mRNAseq.find_ORF(mRNA_seq))

        return mRNAseq(seq_dict[f"{gene_name}_CDS"]['gene_ID'],
                      seq_dict[f"{gene_name}_CDS"]['gene_name'],
                      Sequence.check_seq_type(mRNA_seq),
                      len(mRNA_seq),
                      seq_dict[f"{gene_name}_CDS"]['sp_name'],
                      seq_dict[f"{gene_name}_CDS"]['sbsp_name'],
                      AT_content,
                      TraLa_seq)

    # Method used to translate the open reading frame of a given
    # mRNA sequence to its corresponding amino acid sequence using
    # amino_acid_codons dictionary
    def translate_sequence(ORF):

        mRNAseq.upload_Codons('codon_table.txt')
        codons_dict = mRNAseq.__amino_acid_codons
        codons = ORF
        TraLa_seq = ""



        j = 0
        while (j != len(codons)) and (codons[j] not in codons_dict['O']):
            amino_acid = [key for key, value in codons_dict.items() if codons[j] in value]
            TraLa_seq += amino_acid[0]
            j += 1

        return TraLa_seq

    # The following methods have been used to find the open reading frame
    # of a given mRNA sequence.

    # Method to return the complementary sequence of a given sequence
    def complement(seq):

        complement = ""

        for base in seq:
            if base == 'A':
                complement += 'T'
            if base == 'T':
                complement += 'A'
            if base == 'G':
                complement += 'C'
            if base == 'C':
                complement += 'G'

        return complement

    # Method to find all the possible reading frames of given sequence
    # in one direction of the strand
    def find_reading_frames(mRNA_seq):

        codons = []

        j = 0
        k = 0
        while j < 3:
            i = j
            # print(f"Reading frame {j + 1}")
            while i < len(mRNA_seq):
                if mRNA_seq[i: i + 3] == 'AUG':
                    k = i
                    reading_frame = []
                    while mRNA_seq[k: k + 3] not in ['UAA', 'UAG', 'UGA']:
                        if k == len(mRNA_seq) - 3:
                            reading_frame = []
                            break
                        reading_frame.append(mRNA_seq[k: k + 3])
                        k += 3

                    if len(reading_frame) != 0:
                        # print(f"seq: {reading_frame}")
                        codons.append(reading_frame)
                i += 3
            j += 1

        return codons

    # Method to find open reading frames
    @staticmethod
    def find_ORF(mRNA_seq):

        ORF = []

        # Creating all possible reading frames for,
        # 1) Plus strand
        codons = mRNAseq.find_reading_frames(mRNA_seq)

        for list in codons:

            if len(ORF) == 0:
                for codon in list:
                    ORF.append(codon)

            elif len(list) > len(ORF):
                ORF.clear()
                for codon in list:
                    ORF.append(codon)

        # 2) Minus strand
        # Taking the reverse complement
        complement = mRNAseq.complement(mRNA_seq)
        reverse_complement = complement[::-1]
        codons = mRNAseq.find_reading_frames(reverse_complement)

        for list in codons:

            if len(ORF) == 0:
                for codon in list:
                    ORF.append(codon)

            elif len(list) > len(ORF):
                ORF.clear()
                for codon in list:
                    ORF.append(codon)

        return ORF

class proteinseq(Sequence):

    def __init__(self, UniProt_ID, Reviewed_status, aa_composition, Hydrophobicity, seq):

        self.UniProt_ID = UniProt_ID
        self.Reviewed_status = Reviewed_status
        self.seq_type = Sequence.check_seq_type(seq)
        self.aa_composition = aa_composition
        self.Hydrophobicity = Hydrophobicity

        self.seq_count += 1



    # Method to find the amino acid composition in the given protein
    # sequence
    # Q2 IV
    @staticmethod
    def amino_acid_composition(protein_seq):

        amino_acids = {}

        for amino_acid in protein_seq:
            if amino_acid in amino_acids:
                amino_acids[amino_acid] += 1
            else:
                amino_acids.update({amino_acid : 1})

        return amino_acids

    # Method to calculate the hydrophobicity of a given protein
    def get_Hydrophobicity(protein_seq):

        A = 0
        I = 0
        L = 0
        M = 0
        F = 0
        W = 0
        Y = 0
        V = 0

        for amino_acid in protein_seq:
            if 'A' in amino_acid:
                A += 1
            if 'I' in amino_acid:
                I += 1
            if 'L' in amino_acid:
                L += 1
            if 'M' in amino_acid:
                M += 1
            if 'F' in amino_acid:
                F += 1
            if 'W' in amino_acid:
                W += 1
            if 'Y' in amino_acid:
                Y += 1
            if 'V' in amino_acid:
                V += 1

        return f"{((A + I + L + M + F + W + Y + V)/len(protein_seq))*100} %"

    # The following method automates the creation of a protein sequence
    # object rather than having to manually enter the sequence each
    # time. It uses the sequence information stored in seq_dict
    # dictionary and the relevant gene name and finds the information
    # related to that gene in the dictionary to create the sequence
    # object.
    @staticmethod
    def create_proteinseq_object(seq_dict, protein_name):

        protein_seq = seq_dict[f'{protein_name}_P']['seq']
        aa_composition = proteinseq.amino_acid_composition(protein_seq)
        Hydrophobicity = proteinseq.get_Hydrophobicity(protein_seq)

        return proteinseq(seq_dict[f'{protein_name}_P']['UniProt_ID'],
                          seq_dict[f'{protein_name}_P']['reviewed_status'],
                          aa_composition,
                          Hydrophobicity,
                          protein_seq)


if __name__ == '__main__':

    seq_dict = Sequence.fasta_split('OSDREB_sequences_1.fasta')

    # Q2 I
    DREB1A = DNAseq.create_DNAseq_object(seq_dict, 'DREB1A')

    print(f"Gene_ID: {DREB1A.gene_ID} \n"
          f"Sequence length: {DREB1A.seq_length} \n"
          f"Sequence type: {DREB1A.seq_type}\n"
          f"AT content: {DREB1A.AT_content}")

    # Q2 II & III
    DREB2B_DNA = DNAseq.create_DNAseq_object(seq_dict, 'DREB2B')
    DREB2B_mRNA = mRNAseq.create_mRNAseq_object(seq_dict, 'DREB2B', DREB2B_DNA.TraSc_seq)

    print(f"Sequence length: {DREB2B_mRNA.seq_length}\n"
          f"Sequence type: {DREB2B_mRNA.seq_type}\n"
          f"AT content: {DREB2B_mRNA.AT_content}\n"
          f"mRNA Sequence: \n"
          f"{DREB2B_DNA.TraSc_seq}\n"
          f"Translated Sequence: \n"
          f"{DREB2B_mRNA.TraLa_seq}\n"
          f"Protein sequence length: {len(DREB2B_mRNA.TraLa_seq)}")

    #Q2 IV
    DREB2A_P = proteinseq.create_proteinseq_object(seq_dict, 'DREB2A')
    print(f"UniProt ID: {DREB2A_P.UniProt_ID} \n"
          f"Reviewed status: {DREB2A_P.Reviewed_status} \n"
          f"Sequence type: {DREB2A_P.seq_type}\n"
          f"Amino acid composition: \n{DREB2A_P.aa_composition}\n"
          f"Hydrophobicity: {DREB2A_P.Hydrophobicity}")

    #Q2 V
    print(f"Number of sequences created using the sequence class variable:\n"
          f"{Sequence.seq_count}")





