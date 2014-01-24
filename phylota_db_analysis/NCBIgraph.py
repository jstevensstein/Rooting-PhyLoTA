from os.path import join
import graph_tool.all as gt
from collections import defaultdict, Counter
from itertools import ifilter, combinations, cycle, chain
from pprint import pprint

basepath = '/home/jstevensstein/Desktop/ncbi_db'

node_fields = ["taxid", "parent_taxid", "rank", "embl_code", "division_id",
               "inherited_div_flag", "genetic_code_id", "inherited_gc_flag",
               "mt_gc_id", "inherited_mgc_flag", "genbank_hidden_flag",
               "hidden_subtree_root_flag", "comments"]
name_fields = ["taxid", "name", "unique_name", "name_class"]
keyword_not_clade = ["unclassified", "other", "viral", "viroids", "viruses",
                    "artificial", "unidentified", "endophyte", "scgc", 
                    "libraries","virus", "sample", "metagenome", 
                    "environmental", "uncultured"]


def process_nodes(f):
    n = {}
    r = {}
    c = defaultdict(list)
    i = 0
    for line in f:
        v = [ x.strip() or None for x in line.split("|") ]
        s = dict(zip(node_fields, v[:-1]))
        for k, v in s.items():
            if k.endswith('_flag'):
                s[k] = bool(int(v))
            else:
                try: s[k] = int(v)
                except: pass
        tid = s['taxid']
        r[tid] = s['rank']
        n[tid] = s
        if tid > 1: c[s['parent_taxid']].append(tid)
        i += 1
        print '%s           \r' % i,
    print
    return n, c, r

def process_names(f):
    seen = set()
    synonyms = defaultdict(list)
    accepted = {}
    i = 0
    for line in f:
        v = [ x.strip() or None for x in line.split("|") ]
        v[0] = int(v[0])
        s = dict(zip(name_fields, v[:-1]))
        name = s['name']; uname = s['unique_name']; taxid = s['taxid']
        s['type'] = 'taxonomic_name'
        s['source'] = 'ncbi'
        if uname or name in seen: s['homonym_flag'] = True
        if s['name_class'] == 'scientific name' and (uname or name) not in seen:
            s['primary'] = True
            accepted[taxid] = s
        else:
            s['primary'] = False
            synonyms[taxid].append(s)
        seen.add(uname or name)
        i += 1
        print '%s           \r' % i,
    print
    return accepted, synonyms

def create_taxonomy_graph():
    '''generates a graph_tool graph of the ncbi_hierarchy, and returns that object as well as 3 dictionaries respectively mapping tax id to vertex, name to vertex, and tax id to taxonomic rank'''
    with open(join(basepath, 'nodes.dmp')) as f:
        nodes, ptid2ctid, tid2rank = process_nodes(f)
    with open(join(basepath, 'names.dmp')) as f:
        accepted, synonyms = process_names(f)

    G = gt.Graph()
    v2name = G.new_vertex_property('string')
    G.vertex_properties['name'] = v2name
    v2rank = G.new_vertex_property('string')
    G.vertex_properties['rank'] = v2rank
    v2tid = G.new_vertex_property('int')
    G.vertex_properties['taxid'] = v2tid
    e2tax = G.new_edge_property('bool')
    G.edge_properties['istaxon'] = e2tax
    v2tax = G.new_vertex_property('bool')
    G.vertex_properties['istaxon'] = v2tax

    tid2v = {}
    name2v = {}

    nnodes = len(nodes)
    viter = G.add_vertex(nnodes)

    i = 0
    for tid, d in nodes.iteritems():
        v = G.vertex(i)
        v2tax[v] = 1
        tid2v[tid] = v
        v2tid[v] = tid
        v2rank[v] = tid2rank[tid]
        try:
            name = accepted[tid]['name']
            name2v[name] = v
            v2name[v] = name
        except KeyError:
            print tid
        i += 1
        print '%s           \r' % (nnodes-i),
    print

    i = 0; n = len(ptid2ctid)
    for tid, child_tids in ptid2ctid.iteritems():
        pv = tid2v[tid]
        for ctid in child_tids:
            cv = tid2v[ctid]
            e = G.add_edge(pv, cv)
            e2tax[e] = 1
        i += 1
        print '%s           \r' % (n-i),
    print

    return G, tid2v, name2v, tid2rank


def save_taxonomy_graph(G, filename):
    ''' saves the graph G as filename in the appropriate basepath'''
    G.save(join(basepath, filename)) #must be graphml(xml), gml, or dot format; graphml preferred

def load_taxonomy_graph(filename):
    '''loads the ncbi taxonomy graph saved as filename from the appropriate basepath'''
    print "loading graph"
    G = gt.load_graph(join(basepath, filename))
    n = G.num_vertices
    print "generating tid2v and name2v dictionaries"
    print str(n) + " entries"
    tid2v = {}
    name2v = {}
    v2tid = G.vertex_properties['taxid']
    v2name = G.vertex_properties['name']
    v2rank = G.vertex_properties['rank'] 
    i = 0
    for v in G.vertices():
        tid = v2tid[v]
        name = v2name[v]
        tid2v[tid] = v
        name2v[name] = v

    return G, tid2v, name2v

def fetch_tid2name_name2tid(G):
    '''returns two dictionaries respectively mapping taxon id to scientific name and vice versa'''
    tid2name = {}
    name2tid = {}
    v2tid = G.vertex_properties['taxid']
    v2name = G.vertex_properties['name']
    for v in G.vertices():
        tid = v2tid[v]
        name = v2name[v]
        tid2name[tid] = name
        name2tid[name] = tid
    return tid2name, name2tid

def fetch_mergednodes_dict():
    '''returns a dictionary mapping the two merged taxon ids, the defunct one being the key'''
    f= open(join(basepath, "merged.dmp"))
    merge ={}
    for line in f:
        v = [x.strip() or None for x in line.split("|")]
        merge[int(v[0])] = int(v[1])
    return merge


def fetch_delnodes():
    '''returns a list of taxon ids that have been deleted from the ncbi database from "delnodes.dmp" in the appropriate basepath'''
    f= open(join(basepath, "delnodes.dmp"))
    del_nodes = []
    for line in f:
        del_nodes.append(int(line.split("|")[0].strip()))
    return del_nodes

def add_propmap_taxa_not_clade(G, keyword_not_clade):
    '''adds and internalizes a boolean vertex property map indicating whether taxa are clades'''
    v2name = G.vertex_properties['name']
    taxa_not_clade = G.new_vertex_property('bool')
    G.vertex_properties['not clade flag'] = taxa_not_clade
    for v in G.vertices():
        name = v2name[v].split()
        taxa_not_clade[v] = False
        for kw in keyword_not_clade:
            if kw in name:
                taxa_not_clade[v] = True
                break
    for v in G.vertices():
        if taxa_not_clade[v]:
            seed_taxa_not_clade(v, taxa_not_clade)


def seed_taxa_not_clade(v, taxa_not_clade):
    '''recursively pans through the children of any taxon that is not a clade and flags the others as not clades in the respective boolean property map'''
    taxa_not_clade[v] = True
    for child in v.out_neighbours():
        seed_taxa_not_clade(child, taxa_not_clade)



def is_descendant(ancestor, descendant, tid2v, v2tid):
    '''returns a boolean describing whether or not the first taxon passed is a descendant of the first; can also be treated as a 'is_ancestor function'''
    v = tid2v[descendant]
    while len(list(v.in_neighbours())) != 0:
        if v2tid[v] == ancestor:
            return True
        for n in v.in_neighbours(): #always of length one, if not root
            v = n
    return False

def rank_depth(rank, tid2v):
    v = tid2v[rank]
    depth = 0
    while len(list(v.in_neighbours())) != 0:
        depth += 1
        for n in v.in_neighbours(): #always of length one, if not root
            v = n
    return depth
