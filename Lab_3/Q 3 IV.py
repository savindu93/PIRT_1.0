import networkx as nx

# Q3 IV
# Calculating the number of connections with only a combined score > 0.7 for DREB1A protein

protein_network = nx.read_edgelist("graph_edges.tsv", nodetype = str, data = [('weight', float)])

# PPI with combined scores > 0.7
protein_network_1 = nx.Graph()

for edge in protein_network.edges(data = True):
    # u = edge[0]
    # v = edge[1]
    attr = edge[-1]

    if attr['weight'] > 0.7:
        print(edge)
        protein_network_1.add_edges_from([edge])

print('\n')
print('-'*20)

# Calculating the degree of protein DREB1A (ERF24) with combined score > 0.7

degrees = protein_network_1.degree()

for node, degree in degrees:
    if node == 'ERF24':
        print(f"Degree of {node}(DREB1A) protein with combined score > 0.7: {degree}")
