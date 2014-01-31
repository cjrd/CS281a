"""
HW1 Variable elimination implementation using only std python libraries (python 2.7X)

from the command line, specify clen param
python hw1_eliminate.py --clen=4
assumes that the leaf node is observed with value 1 and calculates the cpd for p(root|leaf)

Author: Colorado Reed (colorado@berkeley.edu)
4 Feburary 2014
"""
import sys
import random
import string
import argparse

class ChainGraph:
    """
    A simple chain graph that implements variable elimination by traversing the input nodes in reverse order
    """
    def __init__(self):
        """
        Chain graph initialization
        TODO consider allowing nodes in constructor
        """
        self.nodes = {}
        self.evidence = {}

    def add_node(self, bnode):
        """
        Add a binary node instance to the graph
        Arguments:
        - `bnode`: binary node to add to graph
        """
        self.nodes[bnode.id] = bnode

    def get_node(self, node_id):
        """
        Get a node by id
        Arguments:
        - `node_id`: the node id
        """
        return self.nodes[node_id]

    def get_nodes(self):
        """
        return the nodes dictionary
        """
        return self.nodes

    def get_evidence(self):
        """
        return the evidence dictionary
        """
        return self.evidence

    def eliminate(self, tnode_id):
        """
        Performs the variable elimination algorithm for the given node id in a
        simple directed chain graph to obtain: p(node_id|conditioned_nodes),
        where conditioned nodes are specified using the
        input_evidence({node_id: evidence}) method
        Arguments:
        - `node_id`: the id of the node that we want to find the CPD for
        """
        # get the leaf node TODO keep track of this structure in a more elegant way
        no_dep = self.nodes.copy()
        for node_id in self.nodes:
            parent = self.nodes[node_id].get_parent()
            if parent is not None:
                del no_dep[parent.get_id()]
        assert(len(no_dep.keys()) == 1)
        leaf = no_dep.keys()[0]

        # find the evidence node in the chain closest to the desired node
        use_leaf = None
        curnode = self.nodes[leaf]
        while curnode.get_id() != tnode_id:
            if curnode.get_id() in self.evidence:
                use_leaf = curnode
            curnode = curnode.get_parent()

        if use_leaf is None:
            raise Exception("No paths between provided evidence and desired target cpd node")

        # do variable elimination
        # TODO test x->y where y is observed
        # TODO test a->x->y where y is observed and x is desired

        # because we're dealing with binary chains,
        # the "influence propagation" is a simple matrix multiplication
        curnode = use_leaf
        # get the prob row corresponding to the evidence
        prop_vec = [curnode.get_cpt()[self.evidence[curnode.get_id()]][:]]
        # move to the parent of the evidence node
        curnode = curnode.get_parent()
        cpt = curnode.get_cpt()
        while curnode.get_id() != tnode_id:
            curnode = curnode.get_parent()
            # TODO this breaks when reaching the root
            new0 = prop_vec[0][0]*cpt[0][0] + prop_vec[0][1]*cpt[1][0]
            new1 = prop_vec[0][0]*cpt[0][1] + prop_vec[0][1]*cpt[1][1]
            prop_vec[0][0] = new0
            prop_vec[0][1] = new1
            # use the child's cpt
            cpt = curnode.get_cpt()

        # add prior probs
        prop_vec[0][0] *= cpt[0][0]
        prop_vec[0][1] *= cpt[0][1]

        # normalize
        Z = prop_vec[0][0] + prop_vec[0][1]
        prop_vec[0][0] /= Z
        prop_vec[0][1] /= Z

        # TODO does this work for chain nodes?
        # TODO FIXME need to incorporate the the cpt of the base node

        # return the cpt
        return prop_vec

    def add_evidence(self, **kwargs):
        """
        Add evidence to the graph using keyword arguments node_id=evidence.
        Input evidence as None to remove evidence
        Arguments:
        - `**kwargs`: input the evidence as keyword arguments, e.g.
        node_id1=observation1, node_id2=observation2, etc
        """
        for nid in kwargs:
            val = kwargs[nid]
            if val is not None:
                self.evidence[nid] = val
            elif nid in self.evidence:
                del self.evidence[nid]

class BinaryNode:
    """
    A node associated with a binary random variable
    """
    def __init__(self, id=None, parent=None, cpt=None):
        """
        Add evidence to the graph using keyword arguments node_id=evidence.
        Input evidence as None to remove evidence
        Arguments:
        - `id`: id of the binary node -- defaults to 8 random letters
        - `parent`: reference to the parent node
        - `cpt`: conditional probability table, where p(node=0|parent=1) would be accessed as cpt[0][1]
        """
        if id is None:
            # assign a random id if we don't input an id
            id = "".join([random.choice(string.letters) for l in xrange(8)])
        self.id = id
        self.parent = parent
        self.cpt = cpt

    def get_cpt_val(self, row, col):
        return self.cpt[row][col]

    # # # Getters and Setters # # #
    def get_parent(self):
        return self.parent

    def get_id(self):
        return self.id

    def get_cpt(self):
        return self.cpt

    def set_parent(self, parent):
        self.parent = parent

    def set_cpt(self, cpt):
        self.cpt = cpt


def build_simple_chain(clength=4, root_cpt=[[0.5, 0.5]], chain_cpt=[[0.6, 0.2], [0.4,0.8]]):
    """
    build a simple chain graph: x1 -> ... -> xclength
    Arguments:
    - `clength`: number of nodes in the chain (default: 4)
    - `root_cpt`: conditional probability table of the root node [[X,Y]] (default: [[0.5, 0.5]])
    - `chain_cpt`: conditional probability table of the chain (non-root) nodes [[X00,X01], [X10, X11]] (default: [[0.6, 0.2], [0.4, 0.8]])
    """
    node_ordering = ["x" + str(cl+1) for cl in xrange(clength)]
    chain_graph = ChainGraph()
    # parent and cpt for x1
    parent_node = None
    in_cpt = root_cpt
    sys.stdout.write("Adding " + str(clength) + " nodes to chain graph...")
    for node_id in node_ordering:
        node = BinaryNode(id=node_id, parent=parent_node, cpt=in_cpt)
        chain_graph.add_node(node)
        # update for next node
        parent_node = node
        in_cpt = chain_cpt
    sys.stdout.write("done\n")

    return chain_graph

def simple_chain_test():
    clen = 4
    rcpt = [[0.5, 0.5]]
    ccpt = [[0.6, 0.2], [0.4,0.8]]
    cgraph = build_simple_chain(clength=clen, root_cpt=rcpt, chain_cpt=ccpt)

    sys.stdout.write("\nexamining properties of constructed graph...")
    nodes = cgraph.get_nodes()
    assert(len(cgraph.get_nodes().keys()) == clen)
    sys.stdout.write("done\n")

    sys.stdout.write("\nchecking cpts...")
    # check the root
    root = cgraph.get_node("x1")
    for i, row in enumerate(rcpt):
        for j, col in enumerate(row):
            assert(rcpt[i][j] == root.get_cpt_val(i,j))
    # check the chain nodes
    for chain_node_id in ["x" + str(cni) for cni in xrange(2,clen+1)]:
        chain_node = cgraph.get_node(chain_node_id)
        for i, row in enumerate(ccpt):
            for j, col in enumerate(row):
                assert(ccpt[i][j] == chain_node.get_cpt_val(i,j))
    sys.stdout.write("done\n")

    # add evidence
    cgraph.add_evidence(x4=1)

    # verify evidence
    assert(cgraph.get_evidence()["x4"] == 1)

    # compute cpt for x1|evidence
    cpval = cgraph.eliminate("x1")
    # TODO make proper asserts
    print "expect 0.4756: " + str(cpval[0][0])
    print "expect 0.5244: " + str(cpval[0][1])


if __name__=="__main__":
    # perform simple variable elimination on a graph
    parser = argparse.ArgumentParser(description='Process command line options for variable elimination.')
    parser.add_argument('--clen', type=int, default=4, help="length of the binary chain (default: 4)")
    parser.add_argument('--mode', dest='mode', type=str, default="infer", help="run mode {infer, test} (default: infer)")
    # TODO add cpt command line specification
    args = parser.parse_args()

    if args.mode == "test":
        print 'starting tests'
        simple_chain_test()
        print 'tests completed'
    elif args.mode == "infer":
        # TODO allow cpts thru arguments (somehow...)
        cgraph = build_simple_chain(clength=args.clen)

        # assume we're traversing the full chain
        cgraph.add_evidence(**{"x" + str(args.clen): 1})
        cpd = cgraph.eliminate("x1")
        print "cpd given evidence is " + str(cpd)
    else:
        raise Exception("--mode option must be 'infer' or 'test', you provided: " + str(args.mode))
