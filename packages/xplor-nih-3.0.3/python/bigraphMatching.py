"""Tools to perform bipartite graph matching.

The code in this module was originally used as described in:

Bermejo, G. A. and Llinas, M.; J. Am. Chem. Soc. 2008, 130, 3797-3805.

Please, cite the above reference if you find it useful. Coded by Guillermo A.
Bermejo.

"""

import copy
import random



# Just a big number. 
infinity = 1.0e100

class Graph:
    """A general graph.

    A graph connects nodes (vertices) by edges (links).  Each edge can also
    have a weight associated with it.  The constructor call is something like:
        g = Graph({'A': {'B': 1, 'C': 2}})
    This makes a graph with 3 nodes, A, B, and C, with an edge of weight 1 from
    A to B, and an edge of weight 2 from A to C.  You can also do:
        g = Graph({'A': {'B': 1, 'C': 2}}, directed=False)
    This makes an undirected graph, so inverse links are also added.  The graph
    stays undirected; if you add more links with g.connect('B', 'C', 3), then
    inverse link is also added.  You can use g.nodes() to get a list of nodes,
    g.get('A') to get a dict of links out of A, and g.get('A', 'B') to get the
    weight of the link from A to B (or None if link doesn't exist).  'Weights'
    can actually be any object at all, and nodes can be any hashable object.

    Extra functionality:
    g.disconnect('A', 'B') removes link from A to B (and the inverse link if
    graph is undirected).
    g.update_weight('A', 'B', 3) sets weight of A->B link to 3, only if link
    exists (the inverse link's weight is also updated in undirected case).
    g.remove_nodes('A', 'B') removes the specified nodes and their incident
    edges (in our example, the updated g contains C only).

    This class has been adapted from: 
    Russell, S. J.; Norvig, P. Artificial Intelligence: a Modern Approach, 2nd
    ed.; Prentice Hall/Pearson Education: Upper Saddle River, NJ, 2003.
    The original class is covered by the MIT licence:
        http://www.opensource.org/licenses/MIT

    """    
    # Numbered comments:
    
    # (1) takes care of links = {A: {B: 1}, B: {}}. If condition in this line
    # is not tested, in the case of an undirected graph, we will make
    # self.dict = {A:{B:1}, B:{}} instead of self.dict = {A:{B:1}, B:{A:1}}.

    # (2) establish connections here by directly modifying self.dict (don't use
    # connect method).

    # (3) existence test of an A-B link must be: self.get(A, B) is not None;
    # simply using self.get(A, B) doesn't work since the link could exist but
    # with zero weight.

    def __init__(self, links=None, directed=True):
        self.dict = {}
        self.directed = directed
        for a in list(links.keys()):
            if a not in list(self.dict.keys()): self.dict[a] = {} # (1)         
            for (b, weight) in list(links[a].items()):
                if b not in list(self.dict.keys()): self.dict[b] = {}
                self.dict[a][b] = weight                    # (2)
                if not self.directed: self.dict[b][a] = weight

    def get(self, a, b=None):
        """Return a link weight or a dict of {node: weight} entries.

        .get(A,B) returns the weight or None if link doesn't exist;
        .get(A) returns a dict of {node: weight} entries (possibly {}) or None
        if A doesn't exist.

        """
        if b is None:
            return self.dict.get(a) # None if a doesn't exist.
        else:
            return self.dict.get(a, {}).get(b)  # None if "a or b" don't exist.

    def connect(self, a, b, weight=1):
        """Add a link from a to b only if both a and b exist in the graph.

        If the graph is undirected also the inverse link is added.  The weight
        of the link(s) is specified by the weight argument.
        
        """
        if a in self.dict and b in self.dict:
            self.dict[a][b] = weight
            if not self.directed:
                self.dict[b][a] = weight

    def disconnect(self, a, b):
        """Remove link from a to b.

        Also the inverse link is removed if the graph is undirected.

        """
        if self.get(a, b) is not None:  # (3)
            del self.dict[a][b]
            if not self.directed:
                del self.dict[b][a]
            
    def update_weight(self, a, b, weight=1):
        """Update weight of link from a to b (if the link already exists).

        Also the weight of the inverse link is updated if the graph is
        undirected.

        """
        if self.get(a, b) is not None: self.connect(a, b, weight)   # (3)

    def nodes(self):
        """Return a list of nodes in the graph."""
        return list(self.dict.keys())

    def remove_nodes(self, *nodes):
        """Remove specified node(s) and its (their) incident edges."""
        for a in nodes:
            del self.dict[a]    # will complain if a not in graph.
            for b in list(self.dict.keys()):
                if a in self.dict[b]: del self.dict[b][a]
        
def UndirectedGraph(links=None):
    """Return a Graph where every edge (including future ones) goes both ways."""
    return Graph(links=links, directed=False)


class Bigraph(Graph):
    """A Graph that represents an undirected bipartite graph (or bigraph).
    
    The constructor call is something like:

    g = Bigraph({'A': {'C': 2},
                 'B': {'C':1, 'D':2},
                 'dummy': {'E': 0, 'F': 0}})

    This makes a bigraph with bipartition (V1, V2), where V1 = ['A', 'B'] and
    V2 = ['C', 'D', 'E', 'F'].  Note how the keys of the input dictionary are
    put into V1, while nodes (vertices) in their associated values into V2.  A
    'dummy' node key can be used to add nodes into V2 which are not connected to
    any node in V1; the 'dummy' node itself is not kept.
    The partition method returns the (V1, V2) tuple representing the
    bipartition.
    
    """
    def __init__(self, links=None):
        temp = []
        for node in list(links.keys()):
            temp.extend(list(links[node].keys()))
        self.nodes1 = [x for x in list(links.keys()) if x != 'dummy']
        self.nodes2 = list(set(temp))

        if 'dummy' in list(links.keys()):
            nodes = list(links['dummy'].keys())
            for node in nodes:
                links[node] = {}
            del links['dummy']

        Graph.__init__(self, links, directed=False)     
 
    def partition(self):
        """Return a tuple with two node lists representing the bipartition."""
        return (self.nodes1, self.nodes2)

    def remove_nodes(self, *nodes):
        """Remove specified node(s) and the corresponding incident edges."""
        Graph.remove_nodes(self, *nodes)
        for node in nodes:  # Update partition.
            if node in self.nodes1: self.nodes1.remove(node)
            if node in self.nodes2: self.nodes2.remove(node)


def build_Knn(bigraph, minimization=False):
    """Return Knn based on an input bigraph.
    
    Knn is a complete bigraph with n nodes in each partition subset.  Knn is
    build from the input bigraph by, if necessary,  adding ("dummy") nodes to
    the smaller node subset until its cardinality equals that of the bigger
    subset (i.e., n).  The edges of the input bigraph ("original edges") are
    kept in Knn, but new zero-weight edges are added to make the output bigraph
    complete.  If minimization=False [default] the weights of the original edges
    are kept.  If minimization=True the weights of the original edges are
    modified so that for each edge: new_weight = W - old_weight, where W is
    bigger (by 1) than the maximum weight in the input bigraph.  As a result,
    edges with big weight in input bigraph have small weight on output, and
    vice versa.    
    
    """
    noedges = 'Input bigraph has no edges'
    # Determine W constant, greater than the maximum weight in input bigraph.
    if minimization:
        weights = []
        for node1 in bigraph.partition()[0]:
            for node2 in bigraph.partition()[1]:
                if bigraph.get(node1, node2) is not None:
                    weights.append(bigraph.get(node1, node2))        
        if not weights: raise noedges
        W = max(weights) + 1    # Greater than max weight by 1.

    links = {}  # Dictionary for constructing Knn.
    for node1 in bigraph.partition()[0]:
        links[node1] = {}
        for (node2, weight) in list(bigraph.get(node1).items()):
            if minimization:
                weight = W - weight
            links[node1].update({node2: weight})

    # Make total number of nodes even (if necessary).
    npad = abs(len(bigraph.partition()[0]) - len(bigraph.partition()[1]))
    for i in range(npad):
        new_node = 'dummy' + str(i+1)
        if len(bigraph.partition()[1]) > len(bigraph.partition()[0]):
            links[new_node] = {}    # Added nodes for output.partition()[0]...
        else:                       # ...are not linked to anything.
            for node in list(links.keys()):
                links[node].update({new_node: 0})   # Added nodes for ...
                # ...output.partition()[0] are zero-weight-linked to all nodes.

    # Add zero-weight edges to make bigraph complete.
    for node1 in list(links.keys()):
        for node2 in bigraph.partition()[1]:
            if bigraph.get(node1, node2) is None:
                links[node1].update({node2: 0})

    return Bigraph(links)
    

def invert_graph_weights(graph):
    """Update link weights of input graph (return nothing).

    For all ij edges with non-zero weight wij, new_wij = W - wij, where W is
    larger than all the wij's.  Zero-weight edges are not updated. After the
    update, the originally maximum weight becomes the minimum (but not zero),
    the minimum (non-zero) becomes the maximum and so on.  This process can be
    thought as a sort of weight "inversion", hence the function name.
    
    """
    # Determine W constant, greater than the maximum weight in input graph.
    weights = []
    for node1 in graph.nodes():
        for node2 in graph.nodes():
            if graph.get(node1, node2) is not None:
                weights.append(graph.get(node1, node2))
    W = max(weights) + 1    # Greater than max weight by 1.

    # Update weights.
    updates = []
    for node1 in graph.nodes():
        for node2 in graph.nodes():
            weight = graph.get(node1, node2)
            if weight is not None:  # If link exists ...
                if weight:          # If weight is non-zero ...
                    new_weight = W - graph.get(node1, node2) # Invert it.
                else:
                    new_weight = weight     # Leave zero weights unchanged.
                updates.append((node1, node2, new_weight))

    for item in updates:
            graph.update_weight(item[0], item[1], item[2])
    return None


class TreeNode:
    """A node in a tree data structure.

    Contains a pointer to the parent node and to the element this node
    represents.  (Based on the search.Node class.)  A tree can be represented by
    a list of nodes, each of them being a tree leaf.  Recursive inquiries of
    parent attribute starting at a leaf, leads to the root node.
    
    """ 
    def __init__(self, element=None, parent=None):
        """Create a search tree Node, derived from a parent."""
        update(self, element=element, parent=parent, depth=0)
        if parent:
            self.depth = parent.depth + 1

    def path(self):
        """Create a list of nodes from the root to this node."""
        x, result = self, [self]
        while x.parent:
            result.append(x.parent)
            x = x.parent
        return result

    def __repr__(self):
        return "<TreeNode %s>" % self.element

                
def best_assign_Kuhn_Munkres(bigraph):
    """Return optimal assignment in bigraph via the Kuhn-Munkres algorithm.

    It follows the treatment by Bondy, J. and Murty R. (1976) "Graph theory
    with applications", New York, North Holland (Chapter 5, Fig. 5.17).  Given
    an input complete bipartite graph (a Bigraph), return a set of (x, y)
    ordered node pairs where x belongs to graph.partition()[0] and y to
    graph.partition()[1], representing the optimal assignment.
    
    """
    # Notation:
    # bigraph.partition()[0] represents X vertex subset in bigraph.
    # bigraph.partition()[1] represents Y vertex subset in bigraph.
    
    # Initial feasible vertex labeling [Eq. 5.12].
    label = {}
    for x in bigraph.partition()[0]:
        label[x] = max([bigraph.get(x, y) for y in list(bigraph.get(x).keys())])
    for y in bigraph.partition()[1]:
        label[y] = 0

    # Equality subgrapgh corresponding to the initial feasible vertex labeling.
    subgraph = build_equal_subgraph(bigraph, label)

    # Set an initial matching, M, on equality subgraph (just one match). 
    # Edges within M are specified by ordered (x, y) vertex pairs where x
    # belongs to subgraph.partition()[0] (X subset) and y to 
    # subgraph.partition()[1] (Y subset). 
    for x in subgraph.partition()[0]:
        if subgraph.get(x):
            y = list(subgraph.get(x).keys())[0] # just one node from x's adjacency.
            M = set([(x, y)])       # Matching M.
            M_saturated_in_X = [x]  # List of M-saturated nodes in X.
            break

    while True:

        # If all vertices in X vertex subset of the equality subgraph are
        # M-saturated, M is an optimal matching [Theorem 5.5]. Stop.
        if len(M_saturated_in_X) == len(subgraph.partition()[0]): return M

##      ### Prints for debuging purposes.
##        print ''
##        print 'M', M
##        print 'Equality subgraph'
##        for x in subgraph.partition()[0]:
##            print x, subgraph.get(x)
##        print ''
##      #################################

        S = set()
        NS = set()
        T = set()

        for x in subgraph.partition()[0]:
            if x not in M_saturated_in_X:
                u = x   # An M-unsaturated vertex in X.
                S.add(u)
                NS.update(list(subgraph.get(u).keys()))
                tree = [TreeNode(u)]# Initialize M-alternating tree rooted at u.
                break

##      ### Prints for debuging purposes.
##        print 'u', u
##        print 'S', S
##        print 'NS', NS
##        print 'T', T
##        print 'NS\T', NS.difference(T)
##      #################################

        while True:

            if NS == T:
                update_labels(bigraph, label, S, T)
                subgraph = build_equal_subgraph(bigraph, label)
                NS = set()  # Reset NS under new subgraph connectivity.
                for x in S:     
                    NS.update(list(subgraph.get(x).keys()))

##              ### Prints for debuging purposes.
##                print '************* NS == T ****************'
##                print 'Equality subgraph'
##                for x in subgraph.partition()[0]:
##                    print x, subgraph.get(x)
##                print ''
##                print 'S', S
##                print 'NS', NS
##                print 'T', T
##                print 'NS\T', NS.difference(T)
##              #################################
                                                            
            y = NS.difference(T).pop()

##          ### Prints for debuging purposes.
##            print 'y', y
##          #################################

            for edge in M:
                if y in edge:   # If y is M-saturated ...
                    z = edge[0]
                    S.add(z)
                    NS.update(list(subgraph.get(z).keys()))
                    T.add(y)
                    tree = new_tree(subgraph, tree, y, z) # Add y and z to tree.
##                  ### Prints for debuging purposes.
##                    print 'z', z
##                    print 'S', S
##                    print 'NS', NS
##                    print 'T', T
##                    print 'NS\T', NS.difference(T)
##                    for node in tree:
##                        print 'tree branch', node.path()
##                  #################################
                    break   
            else:               # ... else if y is M-unsaturated ...
                tree = new_tree(subgraph, tree, y)  # Add only y to tree.
##              ### Prints for debuging purposes.
##                for node in tree:
##                    print 'tree branch', node.path()
##              #################################
                M = M.symmetric_difference(get_path_edges(tree, y))
                M_saturated_in_X = [edge[0] for edge in M]
                break


def build_equal_subgraph(bigraph, label):
    """Return the equality subgraph of bigraph given the vertex labeling."""
    links = {'dummy': {}}
    y_linked = []   # List connected nodes in second partition.
    for x in bigraph.partition()[0]:
        links[x] = {}
        for y in list(bigraph.get(x).keys()):
            weight = bigraph.get(x, y)
            label_sum = label[x] + label[y]
            if safe_eq(label_sum, weight):
                links[x].update({y: weight})

    for y in bigraph.partition()[1]:
        if y not in y_linked:
            links['dummy'].update({y: 0})
                
    return Bigraph(links)

def update_labels(bigraph, label, S, T):
    """Update vertex labeling in input bigraph (no return value)."""
    temp = []
    for x in S:
        for y in bigraph.partition()[1]:
            if y not in T:
                temp.append(label[x] + label[y] - bigraph.get(x, y))
    sigma_label = min(temp)
        
    for node in bigraph.nodes():
        if node in S: label[node] = label[node] - sigma_label
        if node in T: label[node] = label[node] + sigma_label
    return None

def new_tree(bigraph, tree, y, x=None):
    """Return new M-alternating tree by adding elements y and x to input tree.

    Input:  a tree (list of TreeNode inst.) with at list one node (the root).
            y and (optional) x, where root and x belong to
            bigraph.partition()[0] and y to bigraph.partition()[1] of input
            bigraph.  If provided, x is a match of y (i.e., they are paired in
            current matching in bigraph.
    Output: a new tree with added y and (optional) x.

    Rules for tree growing: y is added anywhere in the tree, as long as the
    corresponding edge exists in bigraph.  If provided, x will be a child of y
    in the new tree.

    """
    for leaf in tree:   # Start from every tree leaf.
        node = leaf         # Initialize node.
        new_leaf = None     # Initialize new_leaf.
        while node:     # Descend into branch from current leaf.
            if bigraph.get(node.element, y) is not None: # Attach point found...
                new_leaf = TreeNode(y, node) # ...at current node.
                break
            node = node.parent  # Deeper node in current branch.

        if new_leaf:    # If new leaf is to be added to current tree branch...
            if node == leaf:    # If new leaf doesn't start a new branch...
                tree.remove(leaf)
            if x:
                tree.append(TreeNode(x, new_leaf))  # x always has y as parent.
            else:
                tree.append(new_leaf)

            return tree

def get_path_edges(tree, element):
    """Return a set with M-augmenting path egdes.
    
    The M-augmenting path is the (root-element)-path in input tree.
    
    """
    for leaf in tree:
        if leaf.element == element: path = leaf.path()
    result = set()
    for node1 in path:
        for node2 in path:
            if not node1.depth % 2: 
                if abs(node1.depth - node2.depth) == 1:
                    result.add((node1.element, node2.element))
    return result

    
def second_best_assign(bigraph, M_best):    
    """Return 2nd best assignment on bigraph given the best one.

    An empty set is returned if the second best assignment doesn't exist.  It
    uses the general algorithm outlined in Chegireddy, C.R. & Hamacher, H. W.
    "Algorithms for finding k-best perfect matchings." Discrete Appl. Math. 18:
    155-165, 1987.
    The restriction on this function to be used only for the assignment problem
    (in a bigraph) and not for the more general perfect matching problem (not
    in a bigraph) lies in its reliance on the Kuhn-Munkres algorithm for
    finding best matchings (assignments).
    
    """
    # Steps corresponding to the general algorithm outline in referece paper 
    # are indicated as comments.
    M_best = list(M_best)   # supplied best matching in graph.
    matchings = []
    for i in range(len(M_best)):    # Step 1.
        
        edges = set(M_best[:i]) # Edge set {ej | j < i} in ref. # Step 2.
        new_graph = copy.deepcopy(bigraph)
        for (x, y) in edges:
            new_graph.remove_nodes(x, y)
        new_graph.update_weight(M_best[i][0], M_best[i][1], -infinity) # Step 3.
        
        Ni = best_assign_Kuhn_Munkres(new_graph)  # Step 4.
        weight = edge_weights_sum(new_graph, Ni)    # Can be -infinity.
        
        Ni.update(edges)    # Step 5.
        weight = weight + edge_weights_sum(bigraph, edges)
        
        matchings.append((weight, Ni))

    matchings.sort()    # Sort Ni matchings (i=1,...,len(M_best)) by weight.
    matchings.reverse()
    
    if matchings[0][0] < 0: # If weight of matching is -infinity ...
        return set()
    else:
        return set(matchings[0][1])       
        
        
def k_best_assign(bigraph, K):
    """Return a sorted list with the K best assignments in bigraph.

    The matchings within the output list are sorted in non-increasing order of
    their weights.  It uses the general algorithm outlined by Asratian, A. S.
    et al. (1998) "Bipartite graphs and their applications", Cambridge, U.K.
    New York, Cambridge University Press.
    The restriction on this function to be used only for the assignment problem
    (in a bigraph) and not for the more general perfect matching problem (not
    in a bigraph) lies in its reliance on the Kuhn-Munkres algorithm for
    finding best matchings (assignments).
    
    """
    # References:
    # 1. Asratian, A. et al. (1998) "Bipartite graphs and their applications",
    #    Cambridge, U.K.' New York, Cambridge University Press.
    # 2. Hamacher, H.W. and Queyranne, M. (1986).  "K best solutions to
    #    combinatorial optimization problems." Annals of Operations Research 4:
    #    123-143.
    # 3. Chegireddy, C. R. and Hamacher, H. W. (1987).  "Algorithms for Finding
    #    K-Best Perfect Matchings." Discrete Applied Mathematics 18(2): 155-165.

    # WARNING: this function might crash in cases of degenerate matchings.
    # It tries to solve the problem (cheaply and dirtily) by randomization (see
    # below).  If the function crashes, just run it again until it doesn't.

    M = [set()] * K      # List to hold the K-best matchings (output).
    M[0] = best_assign_Kuhn_Munkres(bigraph)

    if K == 1:
        return [M[0]]   # Pack it in a list to be consistent with the output
                        # when K > 1.

    N = [set()] * K 
    N[0] = second_best_assign(bigraph, M[0])

    I = [set()] * K 
    O = [set()] * K 

    for i in range(1, K):

        max_weight = max([edge_weights_sum(bigraph, an_N) for an_N in N[:i]
                          if an_N])
        for an_N in N[:i]:
            if an_N and edge_weights_sum(bigraph, an_N) == max_weight:
                p = N.index(an_N)

        # Dealing with matching degeneracy:
        # If N[p] already exists the current i+1 best matchings, recalculate it
        # to try to get a different matching with the same weight (degenerate).
        # This bit was introduced by G.B.  I'm sure it can be done more
        # cleverly, e.g., using sensitivity analysis instead of calculating
        # the best matching from scratch using Kuhn-Munkres.
        while N[p] in M:
            N[p] = second_best_restricted(bigraph, I[p], O[p])

        M[i] = N[p].copy()
            

        if M[p].issubset(N[p]):     # This condition appears in ref. 2.
            e = N[p].difference(M[p]).pop() 
            I[i] = I[p].copy()
            I[p].add(e)
            O[i] = O[p].copy()
            O[p].add(e)
        else:
            e = M[p].difference(N[p]).pop()
            I[i] = I[p].copy()
            O[i] = O[p].copy()
            O[i].add(e)
            I[p].add(e)

##        print 'i =', i
##        print 'p =', p
##        print 'e =', e
##        print 'Ip =', I[p], 'Op =', O[p]
##        print 'Ii =',I[i], 'Oi =',O[i]

        # Compute second best matching Np in new Omega(Ip,Op).
        N[p] = second_best_restricted(bigraph, I[p], O[p])

        # Compute second best matching Ni in Omega(Ii,Oi).
        N[i] = second_best_restricted(bigraph, I[i], O[i])

    return M


def second_best_restricted(bigraph, I, O):
    """Return 2nd best assignment in bigraph from restricted solution space.
    
    The restriction on the solution space is determined by Omega(I, O).
    The restriction on this function to be used only for the assignment problem
    (in a bigraph) and not for the more general perfect matching problem (not
    in a bigraph) lies in its reliance on the Kuhn-Munkres algorithm for
    finding best matchings (assignments).
    
    """
    # Generate new restricted graph.
    new_graph = copy.deepcopy(bigraph)
    for (x, y) in I:
        new_graph.remove_nodes(x, y)
    for (x, y) in O:
        new_graph.update_weight(x, y, -infinity)
        
    best = list(best_assign_Kuhn_Munkres(new_graph))  # Best assignment.
    random.shuffle(best)    # Randomly shuffle edges (might help in getting ...
                            # ... distinct second best degenerate assignments).
    secondBest = second_best_assign(new_graph, best)   # Second best assignment.
    if secondBest: secondBest.update(I)
    
    return secondBest


def run_k_best_assign(bigraph, K):
    """Return k_best_assign(bigraph, K).

    This function handles any exeption raised by running k_best_assign and runs
    it again until it terminates properly.  k_best_assign may crash due to the
    randomization process introduced to deal with degeneracy.

    """
    # I don't know what causes the crash in k_best_assign, but we can run it
    # again until it works (I know, notvery elegant). GAB
    while True:
        try:
            return k_best_assign(bigraph, K)
        except:
            continue


def edge_weights_sum(graph, edges):
    """Return the sum of weights of the input edges in graph.

    "edges" can be any container suporting membership test.  An edge within
    edges is a sequence of two of the graph vertices.
    
    """
    weights = []
    for edge in edges:
        weights.append(graph.get(edge[0], edge[1]))
    return sum(weights)

#______________________________________________________________________________
# Some useful stuff.


def update(x, **entries):
    """Update a dict; or an object with slots; according to entries.
    >>> update({'a': 1}, a=10, b=20)
    {'a': 10, 'b': 20}
    >>> update(Struct(a=1), a=10, b=20)
    Struct(a=10, b=20)
    """
    if isinstance(x, dict):
        x.update(entries)   
    else:
        x.__dict__.update(entries) 
    return x

#______________________________________________________________________________
# Safe floating point comparisons (http://www.lahey.com/float.htm).

EPSILON = 0.0000005

def safe_eq(number1, number2):
    """Return "safe" number1 == number2.

    Compares the two input numbers taking into account floating-point
    inexactness.  It can also be used for integers.  Based on "The Perils of
    Floating Point" by Bruce M. Bush (http://www.lahey.com/float.htm).
    
    """
    return abs(number1 - number2) <= abs(number1 + number2) * EPSILON





