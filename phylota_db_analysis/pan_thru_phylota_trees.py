'''run as executable to process the pb.dmp.nr.trees.184 trees file and save them to save_path''' 

import ivy
import graph_tool.all as gt
import NCBIgraph
import tree_to_graph
import convex_colored_subtrees
import string



ncbiG, ncbiG_tid2v, ncbiG_name2v = NCBIgraph.load_taxonomy_graph("db.xml")
tid2name, name2tid = NCBIgraph.fetch_tid2name_name2tid(ncbiG)
ncbiG_v2tid = ncbiG.vertex_properties['taxid']
Eukaryota_tid = name2tid['Eukaryota']
merged_nodes = NCBIgraph.fetch_mergednodes_dict()
del_nodes = NCBIgraph.fetch_delnodes()

ivy.newick.add_label_chars('/#&-')

save_path = '/home/jstevensstein/Desktop/convex_subtrees_3/'

f = open("pb.dmp.nr.trees.184")
data = open(save_path + "tree_data.txt", "w")

i = 0

for line in f:
    try:
        i += 1
        print i
        print "new tree"
        s = line.split('\t')
        tree_name = s[0]
        data.write(tree_name + '\n')
        print "loading tree"
        t = ivy.tree.read(s[1][:-1])
        num_leaves = len(list(t.leaves()))
        data.write("n leaves: " + str(num_leaves) + '\n')
        print "parent tree # leaves = " + str(num_leaves)
        print "converting tree to graph"
        treeG = tree_to_graph.tree_2_G(t)
        tree_to_graph.ammend_treeG_taxid_map(treeG, merged_nodes, del_nodes)
        tree_to_graph.treeG_add_istaxonclade_map(treeG, ncbiG, ncbiG_tid2v)
        tids = tree_to_graph.treeG_to_tids(treeG)
        full_classes = tree_to_graph.tids_to_fullclasses(tids, ncbiG_tid2v, ncbiG_v2tid)
        mrcr = tree_to_graph.fullclasses_to_mrcr(full_classes) #change to ancestor, note that is a taxonomic 'ancestor'
        data.write("mrcr taxid: " + str(mrcr) + '\n')
        if not NCBIgraph.is_descendant(Eukaryota_tid, mrcr, ncbiG_tid2v, ncbiG_v2tid):
            print "not within Eukaryota"
            data.write("not Eukaryota\n")
            continue
        treeG_v2tid = treeG.vertex_properties['taxid']
        tree_to_graph.treeG_add_subclassification_map(treeG, mrcr, ncbiG_tid2v, ncbiG_v2tid)
        tree_to_graph.treeG_add_istaxonclade_map(treeG, ncbiG, ncbiG_tid2v)
        n_leaves, n_unique_tids, n_unclassified_leaves, n_unique_unclassified_tids = tree_to_graph.leaf_stats(treeG, ncbiG_v2tid)
        data.write("n unique taxons: " + str(n_unique_tids) + '\n')
        data.write("n unclassified leaves: " + str(n_unclassified_leaves) + '\n')
        data.write("n unique unclassified taxons: " + str(n_unique_unclassified_tids) + '\n')
        print "finding convex subtrees"
        convex_subtrees, tids_checked = convex_colored_subtrees.extract_convex_subtrees(treeG, merged_nodes, del_nodes, ncbiG_tid2v, ncbiG_v2tid)
        data.write('taxids checked\n')
        data.write(str(tids_checked) + '\n')
        keys = convex_subtrees.keys()
        print "writing extracted trees to file"
        if len(keys) > 0:
            for k in keys:
                file_name = tree_name + '_convextaxid_' + str(k)
                data.write("subtree: " + file_name + '\n')
                tree_file = open(save_path + file_name, "w")
                subtree_data = convex_subtrees[k][1]
                data.write('depth in ncbi taxonomy: ' + str(NCBIgraph.rank_depth(k, ncbiG_tid2v)) + '\n')
                data.write('subtree n leaves: ' + str(subtree_data['n_leaves_tot']) + '\n')
                data.write('subtree taxons:\n')
                data.write(str(subtree_data['tids_and_count']) + '\n')
                data.write('subtree n unclassified leaves: ' + str(subtree_data['n_leaves_not_clade']) + '\n')
                data.write('subtree unique unclassified taxons' + '\n')
                data.write(str(subtree_data['unique_tids_not_clade']) + '\n')
                print "converting to ivy tree and writing"
                newick = ivy.tree.write_newick(convex_subtrees[k][0])
                tree_file.write(newick)
                tree_file.close()
    except:
        print "FAILED"
        data.write('failed')
        #add some catch to describe the error? Note that error usually occurs in the formatting of the tree from pb.dmp.nr.trees.184
    print "resetting"
    data.write('\n')
f.close()