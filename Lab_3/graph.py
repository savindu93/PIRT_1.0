import networkx as nx
sn = nx.Graph() # social network graph

sn.add_edge('Maheshi', 'Dishmi')
sn.add_edge('Maheshi', 'Samadhi')
sn.add_edge('Maheshi', 'Keshari')
sn.add_edge('Sachintha', 'Keshari')

nx.write_gml(sn, "social.gml")

sn_1 = nx.Graph() # social network graph

sn_1.add_edge('Maheshi', 'Dishmi', weight = 0.8)
sn_1.add_edge('Sachintha', 'Dishmi', weight = 0.8)
sn_1.add_edge('Maheshi', 'Samadhi', weight = 0.6)
sn_1.add_edge('Maheshi', 'Keshari', weight = 0.5)
sn_1.add_edge('Sachintha', 'Keshari', weight = 0.3)

nx.write_gml(sn_1, "social_1.gml")

print(nx.shortest_path(sn_1, source = 'Dishmi', target = 'Samadhi'))

# phylotree = nx.DiGraph() # social network graph
#
# phylotree.add_edges_from([("Great Apes", "Humans"), ("Great Apes", "Chimps"), ("Great Apes", "Gorillas")])
#
# nx.write_gml(phylotree, "phylotree.gml")