# Q 3 III
import networkx as nx
import matplotlib.pyplot as plt

# Extracting the edges and weights of each from the .tsv file and
# saving it in a seperate file

f = open("string_interactions_short.tsv", 'r')

graph_edges = ""

for line in f:
    if line.startswith('#'):
        continue
    else:

        list = line.strip().split('\t')
        graph_edges += list[0] + " " + list[1] + " " + list[12] + "\n"

print(graph_edges)

f1 = open("graph_edges.tsv", 'w')
f1.write(graph_edges)

# Creating the graph of the PPI from the above created file
protein_network = nx.read_edgelist("graph_edges.tsv", nodetype = str, data = [('weight', float)])


for edge in protein_network.edges(data = True):
    print(edge)

print('\n\n')
print('-'*20)

# Visualizing the graph
plt.figure(figsize = (8, 6))
nx.draw(protein_network, with_labels = True)
plt.show()

# Calculating the degree of protein DREB1A (ERF24)

degrees = protein_network.degree()

for node, degree in degrees:
    if node == 'ERF24':
        print(f"Degree of {node}(DREB1A) protein: {degree}")




