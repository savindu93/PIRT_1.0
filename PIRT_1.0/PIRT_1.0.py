from PRAT.Protein_Analyzer import PRAT
import streamlit as st



with st.sidebar:

    st.markdown('''
    ## :orange[Protein Information Retriever Tool]
    **PIRT** can retrieve protein data from the Swiss-Prot, Prosite and PDB servers given a set of 
    protein IDs or sequences. It can retrieve information such as the sequences, domains, chains 
    and atoms of the proteins specified.
    
    #### 3rd party packages used  
    [Biopython](https://biopython.org/)  
    [Beautiful Soup](https://beautiful-soup-4.readthedocs.io/en/latest/)
    
    #### Image credits:  
    [3D Illustration of Protein](https://3dproteinimaging.com/wp-content/uploads/2020/08/protein-imager-molecular-illustration-1TOX.jpeg)
    
    :orange[**You can find me at,**] \n
    LinkedIn: [Savindu Weerathunga](https://www.linkedin.com/in/savinduweerathunga/)  
    Twitter: [@savindul_w](https://twitter.com/savindul_w)  
    Source Code: [PIRT 1.0 GitHub Repo](https://github.com/savindu93/PIRT_1.0)
    
    
    ''', unsafe_allow_html = True)

col = st.columns([4, 2])

with col[0]:
    st.markdown("""
    # :blue[PIRT 1.0]
    ## Protein Information Retriever Tool 
    
    """)

with col[1]:
    st.image('PIRT_1.0/protein-illustration-1TOX.png')

tab1, tab2, tab3 = st.tabs(['Retrieve from Swiss-Prot',
                            'Retrieve from Prosite',
                            'Retrieve from PDB'])

with tab1:
    # 1/ Retrieving Swiss-Prot records
    st.subheader('Retrieve Protein Swiss-Prot Records')
    st.markdown(f'''
    
    The following tool retrieves data related to a single or multiple proteins from the
    Swiss-Prot server and outputs the **accession, protein name, organism and the sequence**
    of each protein in **seperate fasta files** which can be downloaded.
    
    {"****"}
    
    ''')
    file = st.file_uploader("**Upload the text file with protein IDs** "
                                "(Protein IDs should be in UniProt Accession format)",
                                type = ["txt"])


    if file and st.button('Retrieve Data', key='retrieve_sp_records'):

        with st.spinner(text = "Retrieving Data"):
            zip_file_path = PRAT.retrieve_SP_records(file)
            #status.update(label = "Retrieved Successfully", state = "complete")

        st.markdown(PRAT.file_downloader(zip_file_path), unsafe_allow_html = True)

    else:
        st.write("Please upload the text file containing the protein IDs.")

with tab2:

    # 2/ Retrieve ExPasy Prosite records containing
    #  domain coordinates, domain IDs, domain names and a description
    #  from the Prosite documentation of the relevant domain
    #  using a text/ fasta file as an input

    st.subheader('Retrieve Domain Information from Prosite')
    st.markdown(f'''

    The following tool retrieves data related to a single or multiple proteins from the
    Prosite server and outputs the following information of each protein in seperate text files
    - **Domain ID**
    - **Domain name**
    - **Domain Coordinates**
    - **A short description of the domain**\n
    
    The proteins can be specified using their UniProt IDs in a text file or sequences in a fasta
    file or can be typed in the given text area.
    
    {"****"}
    
    ''')

    file = st.file_uploader("**Upload a text file containing multiple protein IDs \
                            or fasta file with the protein sequence**",
                                type = ["txt","fasta"])

    user_input = st.text_area("**Enter protein ID or Sequence** "
                              "(When entering multiple IDs enter each new ID in a newline or multiple "
                              "sequences enter each new sequence leaving a blank new line)", "", height = 200)

    if file and st.button('Retrieve Data', key='retrieve_pro_records_f') :

        with st.spinner(text = "Retrieving Data"):
            #zip_file_path = PRAT.retrieve_SP_records(file.name)
            #status.update(label = "Retrieved Successfully", state = "complete")

            domains, zip_file_path = PRAT.retrieve_domain_info_f(file)
            st.markdown(PRAT.file_downloader(zip_file_path), unsafe_allow_html=True)

            for protein, domain_info in domains.items():

                with st.expander(f"{protein}"):
                    st.markdown(domain_info)

    elif user_input and st.button('Retrieve Data', key='retrieve_pro_records_t'):

        with st.spinner(text = "Retrieving Data"):
            #zip_file_path = PRAT.retrieve_SP_records(file.name)
            #status.update(label = "Retrieved Successfully", state = "complete")

            domains, zip_file_path = PRAT.retrieve_domain_info_t(user_input)
            st.markdown(PRAT.file_downloader(zip_file_path), unsafe_allow_html=True)

            for protein, domain_info in domains.items():

                with st.expander(f"{protein}"):
                    st.markdown(domain_info)


    else:
        st.write('''Please upload the text/ fasta file containing the protein IDs
        or Enter the relevant protein ID or sequence''')

with tab3:

    # 1/ Retrieving PDB records/ Chain Sequences
    st.subheader('Retrieve PDB records and Chain Sequences')
    st.markdown(f'''

    The following tool retrieves data related to a single or multiple proteins from the
    PDB server and outputs the,
    - **PDB files**
    - **Chains and their sequences**
    - **Protein sequences given in PDB records**
    
    of each specified protein which can be downloaded.

    {"****"}

    ''')
    file = st.file_uploader("**Upload the text file with PDB protein IDs** ",
                                type = ["txt"])


    if file and st.button('Retrieve Data', key='retrieve_pdb_records'):

        with st.spinner(text = "Retrieving Data"):
            zip_file_path = PRAT.pdb_seq_extractor(file)

            st.markdown(PRAT.file_downloader(zip_file_path), unsafe_allow_html = True)

    else:
        st.write("Please upload the text file containing the protein IDs.")


    # 2/ Retrieving Chain Info from PDB files
    # or Retrieving Chain Info from PDB server when
    # multiple PDB IDs are given in text file
    st.subheader('Retrieve Chain Information ')
    st.markdown(f'''

    The following tool retrieves data related to a single or multiple proteins from the
    PDB server or a PDB file uploaded by the user and outputs the,
    - **No. of chains**
    - **Chain IDs**
    - **No. of residues and no. of atoms in each chain**
    - **Atoms of each residue**

    of each specified protein in seperate text files which can be downloaded.

    {"****"}

    ''')
    file = st.file_uploader("**Upload the PDB file or text file with multiple PDB IDs** ",
                                type = ["ent","pdb","txt"])


    if file and st.button('Retrieve Data', key='retrieve_chain_info'):

        if '.ent' in file.name:
            with st.spinner(text = "Retrieving Data"):
                filepath = PRAT.pdb_chain_extractor_single(file)

                st.markdown(PRAT.file_downloader(filepath), unsafe_allow_html = True)

        elif '.txt' in file.name:
            with st.spinner(text="Retrieving Data"):
                zip_file_path = PRAT.pdb_chain_extractor_multi(file)

                st.markdown(PRAT.file_downloader(zip_file_path), unsafe_allow_html=True)


    else:
        st.write("Please upload the text file containing the protein IDs.")


    # 3/ Retrieving Info on atoms and hetero-residues of given protein in a PDB file
    st.subheader('Retrieve Information on Atoms and Hetero-residues ')
    st.markdown(f'''

    The following tool retrieves data related to a single proteins from a PDB file uploaded 
    by the user and outputs information on the,
    - **Atoms**
    - **Atomic coordinates and other relevant information**
    - **Hetero-residues**
    
    of each specified protein in text files which can be downloaded.

    {"****"}

    ''')
    file = st.file_uploader("**Upload the PDB file** ",
                                type = ["ent","pdb"])


    if file and st.button('Retrieve Data', key='retrieve_atom_hetero_info'):

        with st.spinner(text = "Retrieving Data"):
            filepath = PRAT.pdb_atom_extractor(file)
            data = PRAT.pdb_hetero_extractor(file)

            st.markdown(f"Atomic Info: {PRAT.file_downloader(filepath)}", unsafe_allow_html = True)


            if '.txt' not in data :
                
                st.markdown(f":red[{data}]")

            else:
            
                st.markdown(f"Hetero-residue Info: {PRAT.file_downloader(data)}", unsafe_allow_html = True)






    else:
        st.write("Please upload the text file containing the protein IDs.")








