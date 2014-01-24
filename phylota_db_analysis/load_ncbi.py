'''run as executable to load the NCBI hierarchy and associated structures'''

import ivy
import graph_tool.all as gt
import NCBIgraph
import tree_to_graph
import convex_colored_subtrees
ncbiG, ncbiG_tid2v, ncbiG_name2v = NCBIgraph.load_taxonomy_graph("db.xml")
tid2name, name2tid = NCBIgraph.fetch_tid2name_name2tid(ncbiG)
ncbiG_v2tid = ncbiG.vertex_properties['taxid']
merged_nodes = NCBIgraph.fetch_mergednodes_dict()
del_nodes = NCBIgraph.fetch_delnodes()