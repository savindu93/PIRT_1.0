# Constructing dataframes using the relevant .tsv and .txt files
import pandas as pd


file = "DREB2A_one way_interactions.tsv"
file_1 = "AT_stress_proteins.txt"

column_headers = ['TAIR_Locus','ID','Name','Other_IDs','Plant_Species','DB','Molecule','NCBI_Taxon','GO']

DREB2A_tnw = pd.read_csv(file, sep = "\t")
TAIR_p = pd.read_csv(file_1, sep='\t', header=None)

TAIR_p.columns = column_headers

print(TAIR_p)
print(DREB2A_tnw)

# --------------------------------------------------------------------------------------------------

# Finding the set of known and unknown proteins using the constructed dataframes

# Finding the degree of protein DREB2A
import networkx as nx

dreb2a_nw_edges = ""

for index,row in DREB2A_tnw.iterrows():

    dreb2a_nw_edges += row['#node1'] + " " + row['node2'] + "\n"

f2 = open('dreb2a_edges.tsv', 'w')
f2.write(dreb2a_nw_edges)

dreb2a_nw = nx.read_edgelist('dreb2a_edges.tsv', nodetype = str)

dreb2a_degrees = dreb2a_nw.degree()

for node, degree in dreb2a_degrees:

    if node == 'DREB2A':
        print(f"{node} : {degree}")

# 2 sets to store protein id's of stress related proteins
# and proteins from the network file respectively
stress_p = set()
nw_p = set()

for index,row in TAIR_p.iterrows():

    stress_p.add(row['ID'].lower())

for index, row in DREB2A_tnw.iterrows():

    nw_p.add(row['#node1'].lower())
    nw_p.add(row['node2'].lower())

# Converting all the values in the set back to uppercase
stress_p = {value.upper() for value in stress_p}
nw_p = {value.upper() for value in nw_p}

print(f"Proteins found related to stress: \n"
      f"{stress_p}\n")

print(f"\nProteins in the network:\n"
      f"{nw_p}\n")

# Identifying the proteins that are already known to be involved in
# stress and those are that not known from the network file
known_p = stress_p.intersection(nw_p) # Known proteins
unknown_p = nw_p.difference(stress_p) # Unknown proteins

print(f"\nProteins involved in stress tolerance in the network:\n"
      f"{known_p}\n")

print(f"\nUnknown proteins involved in stress: \n"
      f"{unknown_p}")

print(f"\nNo. of proteins with unknown functions: {len(unknown_p)}\n")


# ----------------------------------------------------------------------------



# Finding the majority voting score for each unknown protein

# Constructing the network with known and unknown proteins without any
# edges/ interactions between unknown proteins (proteins wth unknown functions)

graph_edges = ""
for index,row in DREB2A_tnw.iterrows():

    # Excluding any edges found between unknown proteins
    if row['#node1'] in unknown_p and row['node2'] in unknown_p:

        continue

    elif row['#node1'] in known_p or row['#node1'] in unknown_p:

        graph_edges += row['#node1'] + " " + row['node2'] + " " + str(row['combined_score']) + "\n"


print(graph_edges)
f1 = open("graph_edges.tsv", "w")
f1.write(graph_edges)

p_network = nx.read_edgelist('graph_edges.tsv', nodetype = str, data = [('weight', float)])

for edge in p_network.edges(data = True):
    print(edge)

print('\n')

# Calculating the degree of all proteins in the network
degrees = p_network.degree()

# mv_score stores majority voting scores w.r.t each unknown protein
mv_score = {}

for node, degree in degrees:
    if node in unknown_p:
        mv_score[node] = degree

# Sorting the unknown proteins and their majority voting (mv) score in
# the descending order of its mv score
from collections import OrderedDict

sorted_mv_score = OrderedDict(sorted(mv_score.items(), key = lambda x:x[1], reverse = True))

f3 = open('sorted_mv_score.txt', 'w')
for k,v in sorted_mv_score.items():

    print(f"{k} : {v}")
    f3.write(f"{k} : {v}\n")








