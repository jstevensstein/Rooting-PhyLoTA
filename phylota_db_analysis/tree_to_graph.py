import ivy
import graph_tool.all as gt



def tree_2_G(t): #directly generates a graph from a tree
    '''generates a graph representation of an ivy tree structure including several attributes derived from the label; specialized for use in analyzing phylota's pb.dmp.nr.trees.184'''
    treeG = gt.Graph(directed=False)
    v2name = treeG.new_vertex_property('string')
    treeG.vertex_properties['name'] = v2name
    v2nid = treeG.new_vertex_property('int')
    treeG.vertex_properties['node id'] = v2nid
    v2tid = treeG.new_vertex_property('int')
    treeG.vertex_properties['taxid'] = v2tid
    v2gid = treeG.new_vertex_property('int')
    treeG.vertex_properties['gid'] = v2gid
    edge_length = treeG.new_edge_property('double')
    treeG.edge_properties['length'] = edge_length
    is_leaf = treeG.new_vertex_property('bool')
    treeG.vertex_properties['is leaf'] = is_leaf
    nid2v = {}
    for n in t.iternodes():
        v = treeG.add_vertex(1)
        nid = n.id
        v2name[v] = None
        v2nid[v] = nid
        nid2v[nid] = v
        is_leaf[v] = False
        if (n.isleaf):
            is_leaf[v] = True
            s = n.label
            s = s.split('_')
            ti = s[-1][2:]
            gi = s[-2][2:]
            name = ' '.join(s[:-2])
            v2name[v] = name
            v2gid[v] = gi
            v2tid[v] = ti
    for n in t.iternodes():
        pv = nid2v[n.id]
        for child in n.children:
            cv = nid2v[child.id]
            e = treeG.add_edge(pv, cv)
            edge_length[e] = child.length
    return treeG


def treeG_add_subclassification_map(treeG, mrcr, ncbiG_tid2v, ncbiG_v2tid):
    '''adds the property map of a leaf's taxonomic subclassification (within the tree)'''
    subclass = treeG.new_vertex_property("object")
    treeG.vertex_properties['sub classification'] = subclass
    treeG_tid2v = {}
    treeG_v2tid = treeG.vertex_properties['taxid']
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        if is_leaf[v]:
            tid = treeG_v2tid[v]
            subclass[v] = tid_to_subclassification(tid, mrcr, ncbiG_tid2v, ncbiG_v2tid)



def tid_to_subclassification(tid, mrcr, ncbiG_tid2v, ncbiG_v2tid):
    '''given a taxid and the most recent common rank within a tree graph, returns a list taxon ids representing the taxon's subclassification within the tree graph'''
    v = ncbiG_tid2v[tid]
    classification = []
    classification.append(tid) #note that classification is stored in tids, not ranks, because int comparisons are faster than str cmp
    parent = v.in_neighbours()
    while (len(list(parent)) == 1):
        parent = v.in_neighbours() #find way to remove second call?
        for p in parent: #ought to be obly one entry
            v = p
        tid = ncbiG_v2tid[v]
        if (tid == mrcr):
            break
        classification.append(tid)
        parent = v.in_neighbours()
    classification.reverse() #not abs necessary, but orders classification in intuitive order (highest -> lowest)
    return classification



def tree_to_tids(t, merge_keys, merge, del_nodes):
    '''returns a list of the taxon ids of an ivy tree's leaves'''
    tids = []
    leaves = t.leaves()
    for leaf in leaves:
        s = leaf.label
        tid = int(s.split('_')[-1][2:])
        if tid in del_nodes:
            continue
        if tid in merge_keys:
            tid = merge[tid]
        tids.append(tid)
    return tids


def treeG_to_tids(treeG):
    '''returns a list of the ids of any taxons included in the tree graph'''
    v2tid = treeG.vertex_properties['taxid']
    tids = []
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        if is_leaf[v]:
            tid = v2tid[v]
            tids.append(v2tid[v])
    return tids



def tid_to_fullclassification(tid, ncbiG_tid2v, ncbiG_v2tid):
    '''given tid, returns a list representing the full taxonomic classification of the tid passed'''
    v = ncbiG_tid2v[tid]
    full_classification = []
    full_classification.append(tid) #note that classification is stored in tids, not ranks, because int comparisons are faster than str cmp
    parent = v.in_neighbours()
    while (len(list(parent)) == 1):
        parent = v.in_neighbours() #find way to remove second call?
        for p in parent: #ought to be obly one entry
            v = p
        full_classification.append(ncbiG_v2tid[v])
        parent = v.in_neighbours()
    full_classification.reverse() #not abs necessary, but orders classification in intuitive order (highest -> lowest)
    return full_classification
             

def fullclasses_to_mrcr(full_classes):
    '''given a list of full classifications represented as a list of lists of tids, returns the tax id of the most recent common taxonomic rank'''
    # returns most recent common rank from set of tids
    #full_class is a list of lists of the classification of each organism (by tid)
    working_tid = full_classes[0][0]
    for j in range(0, len(full_classes[0])):
        '''works depsite disparities in length of classification
        because returns at end, after setting working tid appropriately'''
        top_tid = full_classes[0][j]
        for i in range(0, len(full_classes)):
            if full_classes[i][j] != top_tid:
                return working_tid
        working_tid = top_tid
    return top_tid

def tids_to_fullclasses(tids, ncbiG_tid2v, ncbiG_v2tid):
    '''given a list of tids, returns a list of the full classification of those tids'''
    full_classes = []
    for tid in tids:
        v = ncbiG_tid2v[tid]
        full_classes.append(tid_to_fullclassification(tid, ncbiG_tid2v, ncbiG_v2tid))
    return full_classes

def ammend_treeG_taxid_map(treeG, merged_nodes, del_nodes):
    '''ammends the taxon id property map of a tree graph based on a list of deleted nodes and a dictionary mapping defunct tax ids to the taxon ids they have been merged with'''
    v2tid = treeG.vertex_properties['taxid']
    to_remove =[]
    merge_keys = merged_nodes.keys()
    for v in treeG.vertices():
        tid = v2tid[v]
        if tid in del_nodes:
            to_remove.append(v)
        elif tid in merge_keys:
            v2tid[v] = merged_nodes[tid]
    for v in to_remove:
        treeG.remove_vertex(v)
        #possibly add removal if not in ncbi.taxids

def treeG_add_istaxonclade_map(treeG, ncbiG, ncbiG_tid2v):
    '''adds a boolean property map to a tree graph indicating whether each node is contained within a taxon not representing a clade'''
    treeG_v2tid = treeG.vertex_properties['taxid']
    treeG_taxa_not_clade = treeG.new_vertex_property('bool')
    treeG.vertex_properties['not clade flag'] = treeG_taxa_not_clade
    ncbiG_taxa_not_clade = ncbiG.vertex_properties['not clade flag']
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        treeG_taxa_not_clade[v] = False
        if is_leaf[v]:
            tid = treeG_v2tid[v]
            ncbiv = ncbiG_tid2v[tid]
            if ncbiG_taxa_not_clade[ncbiv]:
                treeG_taxa_not_clade[v] = True
    cont = True
    while cont:
        cont = False
        for v in treeG.vertices():
            n = v.out_degree()
            if n > 1: #case internal node
                if not treeG_taxa_not_clade[v]:
                    c = 0 #counter for number of unclassified edges
                    for neighbour in v.all_neighbours():
                        if treeG_taxa_not_clade[neighbour] == 1:
                            c += 1
                    if c >= (n - 1):
                        treeG_taxa_not_clade[v] = True
                        cont = True

def leaf_stats(treeG, ncbiG_v2tid):
    ''' returns a number of statitics about the tree graph; respectively, the number of leaves, the number of unique tax ids, the number of leaves flagged as within a taxon not a clade, and the number of unique such taxon ids'''
    treeG_taxa_not_clade = treeG.vertex_properties['not clade flag']
    treeG_v2tid = treeG.vertex_properties['taxid']
    n_leaves = 0
    n_unclassified_leaves = 0
    unique_tids = []
    unique_unclassified_tids = []
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        if is_leaf[v]:
            tid = ncbiG_v2tid[v]
            tid = treeG_v2tid[v]
            n_leaves += 1
            if tid not in unique_tids:
                unique_tids.append(tid)
            if treeG_taxa_not_clade[v]:
                n_unclassified_leaves += 1
                if tid not in unique_unclassified_tids:
                    unique_unclassified_tids.append(tid)
    n_unique_tids = len(unique_tids)
    n_unique_unclassified_tids = len(unique_unclassified_tids)
    return n_leaves, n_unique_tids, n_unclassified_leaves, n_unique_unclassified_tids