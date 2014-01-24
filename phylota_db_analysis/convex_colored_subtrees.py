import ivy
import graph_tool.all as gt
import tree_to_graph

def extract_convex_subtrees(treeG, merged_nodes, del_nodes, ncbiG_tid2v, ncbiG_v2tid):
    '''returns a dictionary mapping taxon ids to any extracted convex subtrees from the given tree graph'''
    tids_checked = []
    extracted_trees = {}
    taxid = treeG.vertex_properties['taxid']
    subclass = treeG.vertex_properties['sub classification']
    v_in_color = treeG.new_vertex_property('bool') #values should originally be None, after panning through, any remaining None and not False are assumed to be True
    treeG.vertex_properties['v colored'] = v_in_color
    v_extracted = treeG.new_vertex_property('bool')
    treeG.vertex_properties['v extracted'] = v_extracted
    for v in treeG.vertices():
        v_extracted[v] = 1 #maskecd only when 0
    treeG.set_vertex_filter(v_extracted)
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        if is_leaf[v]:
            sub = subclass[v]
            for tid in sub:
                if tid in tids_checked:
                    break
                else:
                    tids_checked.append(tid)
                    print 'checkin tid ' + str(tid) +' for convexity'
                    color_for_tid(tid, treeG) #add some comments, meaning convexity, significance, etc
                    is_convex, root, parent = check_convex(treeG)
                    if is_convex:
                        if check_significance(treeG, root, parent):
                            print 'convex'
                            t = extract_single_convex_subtree(treeG, root, parent)
                            t_mrcr = tree_to_graph.fullclasses_to_mrcr(tree_to_graph.tids_to_fullclasses(tree_to_graph.tree_to_tids(t[0], merged_nodes.keys(), merged_nodes, del_nodes), ncbiG_tid2v, ncbiG_v2tid))
                            extracted_trees[t_mrcr] = t
                            v_extracted[root] = 1
                        print "not significant"
                    else:
                        print "not convex"
    return extracted_trees, tids_checked


def color_for_tid (tid, treeG): #check if a single tid produces a convex colored tree
    #in all cases, within color indicated by red, not in color indicated by black
    '''colors a tree graph appropriately for a given taxon id'''
    v_in_color = treeG.vertex_properties['v colored']
    v_extracted = treeG.vertex_properties['v extracted']
    subclass = treeG.vertex_properties["sub classification"]
    taxa_not_clade = treeG.vertex_properties['not clade flag']
    is_leaf = treeG.vertex_properties['is leaf']
    for v in treeG.vertices():
        v_in_color[v] = 0
        if taxa_not_clade[v]:
            v_in_color[v] = 1  
        if is_leaf[v]:
            sub = subclass[v]
            if tid in sub:
                v_in_color[v] = 1
    cont = True
    while (cont == True):
        cont = recolor_iteration(treeG, v_in_color) #only continue if the recoloring iteration made a change to the coloration, otherwise recoloring is finished


def recolor_iteration(treeG, v_in_color):
    '''represents one iteration through a tree graph to recolor, coloring any uncolored vertex if it has all but one neighbour colored'''
    made_change = False
    for v in treeG.vertices():
        n = v.out_degree()
        if n > 1: #case internal node
            if not v_in_color[v]:
                c = 0 #counter for number of colored edges
                for neighbour in v.all_neighbours():
                    if v_in_color[neighbour] == 1:
                        c += 1
                if c >= (n - 1):
                    v_in_color[v] = 1
                    made_change = True
    return made_change
    '''
    for e in treeG.edges():
        if e_in_color[e] == 1: #only analyze those not marked as false already, true was original setting
            if v_in_color[e.target()] == 0:
                e_in_color[e] = 0
                made_change = True
            if v_in_color[e.source()] == 0:
                e_in_color[e] = 0
                made_change = True
    for v in treeG.vertices():
        if v_in_color[v] == 1 and subclass[v] == None: #only analyze those not marked as false already, true was original setting; don't reanalyze leaves
            e_tot = len(list(v.all_edges()))
            e_not_colored = 0
            for e in v.all_edges():
                if e_in_color[e] == 0:
                    e_not_colored += 1
            if e_not_colored >= e_tot - 1:
                v_in_color[v] = 0
                made_change = True
    '''

def check_significance(treeG, root, parent): #a rooted tree is only meaningful if it has three or more leaves
    '''returns a boolean indicating if the tree graph has three or more leaves, and therefore is significant'''
    taxa_not_clade = treeG.vertex_properties['not clade flag']
    v_in_color = treeG.vertex_properties['v colored']
    is_leaf = treeG.vertex_properties['is leaf']
    n = 0
    for v in treeG.vertices():
        if is_leaf[v]:
            if v_in_color[v]:
                if not taxa_not_clade[v]:
                    n += 1
                    if n > 2:
                        return True
    return False

def decolor_branch(treeG, root, parent, v_in_color):
    is_leaf = treeG.vertex_properties['is leaf']
    '''decolors entire branch up to parent'''
    if is_leaf[root]:
        v_in_color[root] = 0
    else:
        for c in root.out_neighbours():
            if c != parent:
                decolor_branch(treeG, c, root, v_in_color)



def check_convex(treeG):
    '''check the convexity of a colored tree, considering nodes representing taxa that are not clades as both within and without a color (i.e. can be included in a tree, but cannot exclusively make a tree not convex'''
    v_in_color = treeG.vertex_properties['v colored']
    taxa_not_clade = treeG.vertex_properties['not clade flag']
    subclass = treeG.vertex_properties["sub classification"]
    convex_counter = 0
    for v in treeG.vertices():
        if v_in_color[v]:
            if not taxa_not_clade[v]:
                for n in v.out_neighbours():
                    if v_in_color[n] == 0:
                        root = v
                        parent = n
                        convex_counter += 1
                        if convex_counter > 1:
                            return False, None, None
    if convex_counter != 1:
        return False, None, None
    if root.out_degree() == 2: #compensating for case where root is an extraneous node merely connecting two others (rather than three)
        v_in_color[root] = 0
        return check_convex(treeG)
    for n in root.out_neighbours():
        if n != parent:
            if taxa_not_clade[n]: #cannot extract tree where any nodes not-clade are not nested beneath a node that does represent a clade
                decolor_branch(treeG, n, root, v_in_color) #may want to merely decolor the not-clade branch. This is an issue for a node with 4 branches, as 3 are assumed.
                return check_convex(treeG)
    return True, root, parent
    '''
    n = 0 #counter to detect number of colored vertices with uncolored edges. In a convex tree, this is only 1
    root = Nones
    is_convex = False
    for v in treeG.vertices():
        if v_in_color[v] == 1:
            for neighbour in v.all_neighbours(): 
                if not v_in_color[neighbour]:
                        is_convex = True
                        n += 1
                        root = v
                        parent = neighbour
            if n > 1:
                return False, None, None #any vertex with multiple uncolored edges will be uncolored itself, thus this only returns true if the tree is not convex
    #now must pan through until find edge that has no unclassified vertices touching
    true_root = False #boolean representing whether any direct outneighbours are not clades
    while true_root == False and subclass[root] != None:
        true_root = True
        for neighbour in root.out_neighbours():
            if taxa_not_clade[neighbour]:
                true_root = False
                n_true_clade = 0 #must be only one, else cannot continue
                for n in root.out_neighbours():
                    if not taxa_not_clade[n]:
                        n_true_clade += 1
                        if n_true_clade > 1:
                            return False, None, None
                        new_root = n
                parent = root
                root = new_root
    return is_convex, root, parent
    '''


def extract_single_convex_subtree(treeG, root, parent):
    '''Extracts the colored portion of a tree and returns a rooted ivy tree representing that subtree'''
    v_in_color = treeG.vertex_properties['v colored']
    treeG.set_vertex_filter(v_in_color)
    v_extracted = treeG.vertex_properties['v extracted']
    subtree_data = {}
    subtree_data['tids_and_count'] = {}
    subtree_data['unique_tids_not_clade'] = []
    subtree_data['n_leaves_tot'] = 0
    subtree_data['n_leaves_not_clade'] = 0
    t = seed_convex_subtree_from_graph(treeG, root, None, subtree_data)
    t.isroot = True
    v_extracted[root] = 1 #root must be maintained within graph in the case that there are only two convex subtrees comprising the entirety of the graph; need an 'outgroup' as it were
    treeG.set_vertex_filter(v_extracted)
    return [t, subtree_data]

def seed_convex_subtree_from_graph(treeG, curr_vertex, parent_vertex, subtree_data):
    '''recursive function that, given a tree graph, a current vertex and a parent vertex, finds and seeds its children under the new subtree graph'''
    t = ivy.tree.Node()
    v_extracted = treeG.vertex_properties['v extracted']
    subclass = treeG.vertex_properties['sub classification']
    edge_length = treeG.edge_properties['length']
    v_extracted[curr_vertex] = 0 #those with zeros are masked
    is_leaf = treeG.vertex_properties['is leaf']
    if parent_vertex != None:
        t.length = edge_length[treeG.edge(curr_vertex, parent_vertex)]
    if is_leaf[curr_vertex]:
        t.isleaf = True
        subtree_data['n_leaves_tot'] += 1
        taxa_not_clade = treeG.vertex_properties['not clade flag']
        v2name = treeG.vertex_properties['name']
        v2tid = treeG.vertex_properties['taxid']
        v2gid = treeG.vertex_properties['gid']
        name_temp = v2name[curr_vertex].split(' ') #array of words in name
        name_temp = '_'.join(name_temp)
        tid = v2tid[curr_vertex]
        if subtree_data['tids_and_count'].has_key(tid):
            subtree_data['tids_and_count'][tid] += 1
        else:
            subtree_data['tids_and_count'][tid] = 1
        if taxa_not_clade[curr_vertex]:
            subtree_data['n_leaves_not_clade'] += 1
            if tid not in subtree_data['unique_tids_not_clade']:
                subtree_data['unique_tids_not_clade'].append(tid)
        gid = v2gid[curr_vertex]
        label = name_temp + '_' + 'gi'  + str(gid) + '_' + 'ti' + str(tid)
        t.label = label
    else:
        t.isleaf = False
        for child in curr_vertex.out_neighbours():
            if child != parent_vertex:
                t.add_child(seed_convex_subtree_from_graph(treeG, child, curr_vertex, subtree_data))
    return t