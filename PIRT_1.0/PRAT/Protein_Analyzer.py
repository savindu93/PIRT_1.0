from Bio import SwissProt
from Bio import ExPASy
from Bio.ExPASy import ScanProsite
import re
from Bio import SeqIO
import os
import zipfile
import base64
import pandas as pd
from io import StringIO
import streamlit as st
from urllib.parse import urlencode

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)

class PRAT:

    # Method to retrieve a set of SwissProt records with IDs
    # given as a list

    def retrieve_SP_records(file):

        try:
            IDs = file.getvalue().decode("utf-8").split("\n")
            protein_data = []

            for ID in IDs:
                handle = ExPASy.get_sprot_raw(ID.strip())
                protein_data.append(SwissProt.parse(handle))

            filenames = []

            for records in protein_data:
                for record in records:
                    file_name = f'{record.accessions[0]}.fasta'
                    filenames.append(file_name)
                    organism = f"{record.organism.split(' ')[0]}_{record.organism.split(' ')[1]}"
                    with open(file_name,'w') as file:
                        file.write(f">{record.accessions[0]} {record.entry_name} {organism}\n"
                                f"{record.sequence}\n")

            filepaths = [os.path.abspath(filename) for filename in filenames]

            zip_file_path = "fasta_files.zip"
            with zipfile.ZipFile(zip_file_path, "w") as zipf:
                for filepath in filepaths:
                    zipf.write(filepath, os.path.basename(filepath))

            return zip_file_path
        
        except Exception as e:
            st.error(f"An error occurred: {e}. \n"\
            "Recheck your input IDs in the text file and\n"\
            "retry.")

    def file_downloader(filepath):

        with open(filepath, 'rb') as f:
            data = f.read()

        b64 = base64.b64encode(data).decode()
        href = f"<a href = 'data:application/zip;base64,{b64}' download = '{os.path.basename(filepath)}'> Download Files </a>"
        return href

    #----------------------------------------------------------------------------------------------------------------------

    # Method to retrieve ExPasy Prosite records containing
    # domain coordinates, domain IDS, domain names and Prosite
    # documentation of the relevant domain using a text/ fasta file
    # as an input

    # Function to retrieve domain info from the scanned results from Prosite
    def scan_domain_info(result):

        domain_info = ""

        i = 0

        for record in result:

            i += 1

            domain_id = record['signature_ac']
            domain_coords = f"{record['start']} - {record['stop']}"

            # Obtain additional domain information related to the domain
            # from Prosite documentation
            import requests
            response = requests.get("https://prosite.expasy.org/" + domain_id)
            # print(response.text)

            from bs4 import BeautifulSoup
            prosite_DOC_Id = "(PDOC[0-9]{5})"
            pattern = re.compile(prosite_DOC_Id)
            soap = BeautifulSoup(response.content, 'html.parser')
            ps_general_info = soap.find('table', class_="type-1").get_text()
            prosite_ID = pattern.findall(ps_general_info)
            # print(prosite_ID)

            domain_name = soap.find_all('td', attrs={'property': 'schema:description'})
            domain = ''
            for names in domain_name:
                domain = names.text.strip('\n\t').replace('domain profile.', '')

            response = requests.get("https://prosite.expasy.org/" + prosite_ID[0])
            # print(response.text)
            soap = BeautifulSoup(response.content, 'html.parser')
            spans = soap.find_all('span', attrs={'property': 'schema:description'})

            # print(f"Domain: {i + 1}\n" \
            #       f"Domain Name: {domain}\n"
            #       f"Domain ID: {domain_id}\n"
            #       f"Domain Coordinates: {domain_coords}\n")
            #
            # for span in spans:
            #     print(span.text)

            domain_info += f"{'*' * 50}\n" \
                           f"Domain: {i} |\n" \
                           f"Domain Name: {domain} |\n" \
                           f"Domain ID: {domain_id} |\n" \
                           f"Domain Coordinates: {domain_coords} \n"

            for span in spans:
                domain_info += f"{span.text}\n"

        return domain_info

    # Retrieve domain information when the protein IDs/ sequences are given in file
    def retrieve_domain_info_f(file):

        domains = {}

        # Retrieve domain info of proteins given in a text file as a list of UniProt accs
        if '.txt' in file.name:

            prot_accs = file.getvalue().decode("utf-8").split("\n")

            for uniprot_acc in prot_accs:

                # Obtain domain ID, name and coordinates
                # relevant to the given UniProt/Swiss-Prot ID
                # from Prosite

                handle = ScanProsite.scan(seq=uniprot_acc.strip())
                result = ScanProsite.read(handle)
                #print(result)

                if result:
                    domain_info = PRAT.scan_domain_info(result)
                else:
                    domain_info = "No domain hits have been found for this protein within the PROSITE server"

                domains[f'UniProt Acc: {uniprot_acc}'] = domain_info

        # Retrieve domain info of proteins given in a fasta file
        else:
            file_handle = StringIO(file.getvalue().decode("utf-8"))
            records = SeqIO.parse(file_handle, format = 'fasta')

            for record in records:

                sequence = record.seq

                # The regular expression pattern to identify a
                # standard UniProt ID (* cited from UniProt)
                uniprot_id = "([OPQ][0-9][A-Z0-9]{3}[0-9])|([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
                pattern = re.compile(uniprot_id)
                match = pattern.findall(record.id)

                uniprot_acc = match[0][0]

                if match == []:
                    # If the protein UniProt ID is not given
                    # find domain info using sequence

                    handle = ScanProsite.scan(seq = sequence)
                    result = ScanProsite.read(handle)

                    domain_info = PRAT.scan_domain_info(result)

                else:
                    # Obtain domain ID, name and coordinates
                    # relevant to the given UniProt/Swiss-Prot ID
                    # from Prosite

                    handle = ScanProsite.scan(seq = uniprot_acc)
                    result = ScanProsite.read(handle)

                    domain_info = PRAT.scan_domain_info(result)


                if domain_info == '':
                    domain_info += "No domain hits have been found for this protein within the PROSITE server"


                domains[f'UniProt Acc: {uniprot_acc} | Protein Name: {record.name}'] = domain_info

        # Create files to download them together in a zipped file
        filenames = []
        for protein, domain_info in domains.items():

            filename = f"{protein.split('|')[0]}_Domains"
            filenames.append(filename)

            with open(filename, 'w', encoding = 'utf-8') as file:
                file.write(domain_info)

        filepaths = [os.path.abspath(filename) for filename in filenames]

        zip_file_path = "domain_files.zip"
        with zipfile.ZipFile(zip_file_path, "w") as zipf:
            for filepath in filepaths:
                zipf.write(filepath, os.path.basename(filepath))

        return domains, zip_file_path

    # Retrieve domain information when the protein IDs/ sequences are given as text input
    def retrieve_domain_info_t(text):

        domains = {}
        i = 0

        uniprot_id = "([OPQ][0-9][A-Z0-9]{3}[0-9])|([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
        pattern = re.compile(uniprot_id)
        match = pattern.findall(text)

        if match:
            inputs = text.split("\n")
        else:
            inputs = text.split("\n\n")

        for input in inputs:

            i += 1
            # Obtain domain ID, name and coordinates
            # relevant to the given UniProt/Swiss-Prot ID
            # from Prosite

            handle = ScanProsite.scan(seq=input)
            result = ScanProsite.read(handle)

            if result:
                domain_info = PRAT.scan_domain_info(result)

            else:
                domain_info = "No domain hits have been found for this protein within the PROSITE server"

            uniprot_id = "([OPQ][0-9][A-Z0-9]{3}[0-9])|([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
            pattern = re.compile(uniprot_id)
            match = pattern.findall(input)

            if match:
                domains[f'Protein: {match[0][0]}'] = domain_info
            else:
                domains[f'Protein: {i}'] = domain_info

        # Create files to download them together in a zipped file
        filenames = []
        for protein, domain_info in domains.items():
            filename = f"{protein.split('|')[0]}_Domains"
            filenames.append(filename)

            with open(filename, 'w', encoding='utf-8') as file:
                file.write(domain_info)

        filepaths = [os.path.abspath(filename) for filename in filenames]

        zip_file_path = "domain_files.zip"
        with zipfile.ZipFile(zip_file_path, "w") as zipf:
            for filepath in filepaths:
                zipf.write(filepath, os.path.basename(filepath))

        return domains, zip_file_path


    #----------------------------------------------------------------------------------------------------------------------

    # Method to retrieve PDB files when PDB IDs are given
    # which should also output a single FASTA file with all
    # PDB sequences


    def pdb_seq_extractor(file):

        # Variable to store the pdb sequences
        # in fasta format

        amino_acid_dict = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Gln': 'Q',
                           'Glu': 'E', 'Glx': 'Z', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                           'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y',
                           'Val': 'V'}

        filenames = []
        filepaths = []
        seq = ''

        IDs = file.getvalue().decode("utf-8").split("\n")

        for ID in IDs:

            pdb_id = ID.strip()

            # Download and Read the PDB file
            from Bio.PDB.PDBList import PDBList
            pl = PDBList()
            protein_file = pl.retrieve_pdb_file(pdb_code = pdb_id, file_format = 'pdb')
            filepaths.append(protein_file)


            from Bio.PDB.PDBParser import PDBParser
            parser = PDBParser(PERMISSIVE = 1)
            struct = parser.get_structure(pdb_id, protein_file)

            models = [model for model in struct.get_models()]


            # Extract the residues for the protein
            seq_res = [res.resname for res in struct.get_residues() if res.resname.lower() in [aat.lower() for aat in amino_acid_dict.keys()]\
                        and struct.get_models() == models[0]] # Only amino acid residues of the 1st model are included, hetero residues are removed 
            
            seq += f'\n>{struct.get_full_id()[0]}\n'

            for res in seq_res:

                seq += [aa for aat,aa in amino_acid_dict.items() if aat.lower() == res.lower()][0]

            seq += '\n'

            # Extract residues for each chain in the protein
            chains = [value for value in struct.get_chains()]

            i = 0
            chain_info = ""
            chain_ids = []

            for chain in chains:

                if chain.get_id() in chain_ids:
                    break
                else:
                    chain_ids.append(chain.get_id())

                chain_seq = ""
                chain_res = [res.resname.lower() for res in chain.get_residues() if res.resname.lower() in [aat.lower() for aat in amino_acid_dict.keys()]]
                #print(chain_res)

                for res in chain_res:

                    chain_seq += [aa for aat,aa in amino_acid_dict.items() if aat.lower() == res][0]

                chain_info += f"\n>{chain.get_full_id()[0]}| Chain: {chain.get_id()}\n " \
                                f"{chain_seq}\n"

                # print(f"\n>{chain.get_full_id()[0]}| Chain: {chain.get_id()}\n "
                #       f"{chain_seq}\n")

                i += 1

            # Output the chain sequences for a specific protein into a file
            # (seperate files for each protein)
            with open(f"{struct.get_full_id()[0]}_chain_seq.txt", 'w') as file:

                file.write(chain_info)
                filenames.append(f"{struct.get_full_id()[0]}_chain_seq.txt")


        # Write the pdb sequences into an output fasta file
        # (single file for all pdb sequences)
        with open('pdb_seq.fasta', 'w') as file:

            file.write(seq)
            filenames.append('pdb_seq.fasta')

        # Create filepaths for the above newly formed files
        for filename in filenames:

            filepaths.append(os.path.abspath(filename))

        zip_file_path = "pdb_files.zip"
        with zipfile.ZipFile(zip_file_path, "w") as zipf:
            for filepath in filepaths:
                zipf.write(filepath, os.path.basename(filepath))

        return zip_file_path


    # ----------------------------------------------------------------------------------------------------------------------

    # Method to output the no.of chains, and atoms
    # for each residue of a protein given in a single pdb file

    def pdb_chain_extractor_single(file):

        stringio_file = StringIO(file.getvalue().decode("utf-8"))


        ex = r'pdb(\w{4}).ent'
        pattern = re.compile(ex)

        pdb_id = re.findall(pattern, file.name)[0]
        print(pdb_id)

        # Read the PDB file
        from Bio.PDB.PDBParser import PDBParser
        parser = PDBParser(PERMISSIVE = 1)
        struct = parser.get_structure(pdb_id, stringio_file)

        no_models = len([model for model in struct.get_models()])

        chains = [value for value in struct.get_chains()]

        print(f"Number of chains: {len(chains)}\n")

        # Variable that includes all the data that will be printed out
        # to the external text file

        if no_models > 1:

            data = f"Protein's Chain Information\n" \
                    f"Protein : {pdb_id}\n\n"\
                    f"Number of models: {no_models}\n\n"
        else:

            data = f"Protein's Chain Information\n" \
                    f"Protein : {pdb_id}\n\n"\
                    f"Number of chains: {len(chains)}\n\n"

        # Extract the information related all the protein chains
        i = 0 # counter for number of models

        for chain in chains:

            residues = [value for value in chain.get_residues()]

            no_atoms = 0
            atoms = []
            for residue in residues:
                residue_atoms = [value for value in residue.get_atoms()]
                no_atoms += len(residue_atoms)

                atom_labels = []
                for atom in residue_atoms:
                    atom_labels.append(atom.get_id())

                atoms.append([residue.get_resname(),atom_labels])

            print(f"Number of atoms: {no_atoms}")

            df = pd.DataFrame(atoms, columns = ['Residue','Atoms'])

            # Add the information to the data variable that will be printed out
            # to the text file

            if no_models > 1:

                i += 1
                data += f"Chain: {chain.get_id()}\n" \
                        f"Model: {i}\n"
            
            else:

                data += f"Chain: {chain.get_id()}\n"
                
            data += f"Number of residues: {len(residues)}\n" \
                    f"Number of atoms: {no_atoms}\n\n" \
                    f"{df.to_string()}\n\n"


        # Write the chain info into an output text file
        with open(f"{pdb_id}_Chain_info.txt",'w') as file:

            file.write(data)
            filename = f"{pdb_id}_Chain_info.txt"

        # Create filepaths for the above newly formed files
        filepath = os.path.abspath(filename)


        return filepath

    # ----------------------------------------------------------------------------------------------------------------------

    # Modified Method to output the no.of chains, and atoms
    # for each residue of proteins given in a text file as pdb IDs

    # Function to extract info from individual chains/ models
    def extract_chain_info(chains, model):

        import pandas as pd
        pd.set_option("display.max_rows", None)

        data = ''

        # Extract the information related to all the protein chains/ models
        for chain in chains:

            residues = [value for value in chain.get_residues()]

            print(f"Chain: {chain.get_id()}\n"
                  f"Number of residues: {len(residues)}")

            no_atoms = 0
            atoms = []

            # Extracting the atoms and their labels w.r.t a specific
            # residue
            for residue in residues:
                residue_atoms = [value for value in residue.get_atoms()]
                no_atoms += len(residue_atoms)

                atom_labels = []
                for atom in residue_atoms:
                    atom_labels.append(atom.get_id())

                atoms.append([residue.get_resname(), atom_labels])

            print(f"Number of atoms: {no_atoms}")

            df = pd.DataFrame(atoms, columns=['Residue', 'Atoms'])
            print(df)

            # Add the information to the data variable that will outputted
            # to the text file
            data += f"Chain: {chain.get_id()}\n" \
                    f"Number of residues: {len(residues)}\n" \
                    f"Number of atoms: {no_atoms}\n" \
                    f"\n{df.to_string()}\n\n"

            if model == True:
                break

        return data

    def pdb_chain_extractor_multi(file):


        filenames = []
        filepaths = []

        chain_data = {}

        IDs = file.getvalue().decode("utf-8").split("\n")

        for ID in IDs:

            st.write(ID)
            pdb_id = ID.strip()

            # Download and Read the PDB file
            from Bio.PDB.PDBList import PDBList
            pl = PDBList()
            protein_file = pl.retrieve_pdb_file(pdb_code = pdb_id, file_format = 'pdb')

            from Bio.PDB.PDBParser import PDBParser
            parser = PDBParser(PERMISSIVE = 1)
            struct = parser.get_structure(pdb_id, protein_file)

            exp_method = struct.header['structure_method'].split(';')

            # Variable that includes all the data that will be printed out
            # to the external text file
            data = f"Protein information\n" \
                    f"Protein : {pdb_id}\n" \
                    f"Experimental Method: {exp_method[0]}\n\n"

            chains = [value for value in struct.get_chains()]

            # If the protein structure has been determined using any of the 3 exp methods given below then such 
            # records will have more than one model for a particular protein chain. By setting model as true
            # only one record out of the multiple models will be printed out.

            if any(method in exp_method for method in ["solution scattering","infrared spectroscopy","solution nmr"]) :

                models = [value for value in struct.get_models()]
                data += f"Number of models: {len(models)}\n\n"

                model = True

                data += PRAT.extract_chain_info(chains,model)

            else:

                data += f"Number of chains: {len(chains)}\n\n"

                model = False
                data += PRAT.extract_chain_info(chains,model)

            chain_data[f'{pdb_id}'] = data

            # Write the chain info into an output text file
            with open(f"{pdb_id}_Chain_info.txt", 'w') as file:

                file.write(data)
                filenames.append(f"{pdb_id}_Chain_info.txt")

        # Create filepaths for the above newly formed files
        for filename in filenames:
            filepaths.append(os.path.abspath(filename))

        zip_file_path = "Chain_info.zip"
        with zipfile.ZipFile(zip_file_path, "w") as zipf:
            for filepath in filepaths:
                zipf.write(filepath, os.path.basename(filepath))

        return zip_file_path

    # Method to extract the atomic coordinates of each residue of a given protein in PDB file format
    def pdb_atom_extractor(file):

        stringio_file = StringIO(file.getvalue().decode("utf-8"))

        ex = r'pdb(\w{4}).ent'
        pattern = re.compile(ex)

        pdb_id = re.findall(pattern, file.name)[0]

        # Variable to store the information on the atoms and their coordinates
        atom_coords = []


        pdb_file = stringio_file.read().split("\n")
        #st.write(pdb_file)


        for line in pdb_file:

            if line.startswith('ATOM'):

                row_data = line.strip('\n').split()

                pattern = re.compile(r"\d.\d{5}.\d{2}")
                col_data = row_data[9]


                if re.match(pattern, col_data):

                    row_data.append(row_data[10])
                    row_data[9] = ''.join(list(col_data)[0:4])
                    row_data[10] = ''.join(list(col_data)[4:])

                atom_coords.append(row_data[1:12])

        columns = ['Atom No.','Atom Name','Amino Acid','Chain ID','AA No.','X','Y','Z','Occupancy','Temperature Factor','Element Symbol']
        df = pd.DataFrame(atom_coords, columns = columns)
        df.set_index(df.columns[0], inplace = True)


        df = df.to_string()

        with open(f'{pdb_id}_atom_info.txt', 'w') as file:

            data = f"Atomic Information\n" \
                   f"Protein: {pdb_id}\n\n" \
                   f"The following file gives information on the atoms and relevant information such as the orthogonal coordinates\n " \
                   f"in X, Y, and Z axes in Angstroms, the occupancy and temperature factor for each residue of given protein.\n" \
                   f"Any integer following the letter in the Element Symbol column refers to the charge of that atom. \n" \
                   f"AA No.: Amino Acid No.\n\n" \
                   f"{df}"

            file.write(data)
            file_name = f'{pdb_id}_atom_info.txt'


            # Create filepath for the above newly formed file
            filepath = os.path.abspath(file_name)


        return filepath


    # Method to extract the hetero residues of a given protein
    def pdb_hetero_extractor(file):
        
        stringio_file = StringIO(file.getvalue().decode("utf-8"))

        ex = r'pdb(\w{4}).ent'
        pattern = re.compile(ex)

        pdb_id = re.findall(pattern, file.name)[0]
        print(pdb_id)

        # Variable to store the information on the hetero-residues
        hetero_residue = []
        het_info = []

        pdb_file = stringio_file.read().split("\n")


        for line in pdb_file:

            if line.split(" ")[0] == 'HET':
                hetero_residue.append(line.strip('\n').split()[1:])

            elif line.split(" ")[0] == 'HETNAM':

                hetnam = line.strip('\n').split(" ")[1:]
                print(hetnam)

                if len(hetnam) > 2:
                    if hetnam[0] == '2':
                        het_name = " ".join(hetnam[2:])

                        het_info[-1][1:] = []
                        het_info[-1].append(het_name)

                    else:
                        het_name = " ".join(hetnam[1:])

                        hetnam[1:] = []

                        hetnam.append(het_name)
                        het_info.append(hetnam)

                else:
                    het_info.append(hetnam)


        hetero_residue_columns = ['Het ID', 'Chain ID', 'Sequence No.', 'No. of Hetero Residues']
        het_info_columns = ['Het ID', 'Het Name']

        df = pd.DataFrame(hetero_residue, columns = hetero_residue_columns)
        df.set_index(df.columns[0], inplace=True)

        df_1 = pd.DataFrame(het_info, columns = het_info_columns)
        df_1.set_index(df_1.columns[0], inplace=True)

        df_0 = df.to_string()
        df_1 = df_1.to_string()

        if df.empty:

            error = "No hetero-residues were found for the specified protein"

            return error

        else:

            with open(f'{pdb_id}_hetero_info.txt', 'w') as file:

                data = f"Hetero-Residue Information\n" \
                    f"Protein: {pdb_id}\n\n" \
                    f"The following file gives information on the hetero-residues of the protein of interest. Information \n" \
                    f"such as the name, number of each hetero residue, the chain and the location of it in the sequence are \n" \
                    f"given.\n\n" \
                    f"{df_0}\n\n" \
                    f"\n{df_1}"

                file.write(data)
                file_name = f'{pdb_id}_hetero_info.txt'


                # Create filepath for the above newly formed file
                filepath = os.path.abspath(file_name)

            return filepath




    # ----------------------------------------------------------------------------------------------------------------------









