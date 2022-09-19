from math import log10
from typing import List, Dict

def opposite_color(color):
    return 'A' if color == 'R' else 'R'

class Gene:
    __attr__ = "in_repr", \
               "in_act", \
               "out_repr", \
               "out_act", \
               "out_m", \
               "active_out_m", \
               "active_in_m", \
               "in_m", \
               "name", \
               "self_act", \
               "self_repr", \
               "edgelist", \
               "ind", \
               "state", \
               "in_act_neighbors", \
               "in_repr_neighbors", \
               "out_act_neighbors", \
               "out_repr_neighbors", \
               "R_j_estimate"\
               "bconnect"  # dict with binode connections

    # in_act_neighbors, in_repr_neighbors, out_act_neighbors, our_repr_neighbors <- can be estimated from the dataset
    # additional attributes for Gene: in_repr_repr, in_repr_act s.t. in_repr = in_repr_repr + in_repr_act <- needs the node coloring to estimate that

    # method for estimation add.attr:
    # estimate_bi_node_connections(self, node_coloring): \\O(dev(self))

    def estimate_bi_node_connections(self, node_coloring):
        res = {'out_A_A': 0,
               'out_R_A': 0,
               'out_A_R': 0,
               'out_R_R': 0,
               'in_A_A': 0,
               'in_R_A': 0,
               'in_A_R': 0,
               'in_R_R': 0,
               'self_A': 0,
               'self_R': 0
               }  # type of connection [out/in] [Activation edge/Repression edge] [ToActivator/toRepressor] or [fromAct/fromRepr]
        res['self_A'] += self.self_act  # add self edges from the Gene
        res['self_R'] += self.self_repr
        all_neighbors = [self.out_act_neighbors, self.out_repr_neighbors, self.in_act_neighbors, self.in_repr_neighbors]
        keys = ['out_A_', 'out_R_', 'in_A_', 'in_R_']
        for dir_neighbors, key in zip(all_neighbors, keys):
            for u in dir_neighbors:
                color = node_coloring[u.ind]
                res[key + color] += 1
        self.bconnect = res

    def update_connections(self, new_color):
        """Function to update the binode connections of the neighbors of the vertex. Run when the color of the current vertex is changed to the new_color and all of it neighbors need to adjust their corresponding counters bconnect
        """
        all_neighbors = [self.out_act_neighbors, self.out_repr_neighbors, self.in_act_neighbors, self.in_repr_neighbors]
        keys = ['out_A_', 'out_R_', 'in_A_', 'in_R_']
        for dir_neighbors, key in zip(all_neighbors, keys):
            for u in dir_neighbors:
                u.bconnect[key + new_color] += 1
                u.bconnect[key + opposite_color(new_color)] -= 1

    def __str__(self):
        return "Gene's name: {0}, out_act: {1}, out_repr: {2}, in_act: {3}, in_repr: {4}, Total out edges: {5}, " \
               "Total in edges: {6}".format(self.name, self.out_act, self.out_repr, self.in_act, self.in_repr,
                                            self.out_m, self.in_m)

    def print_edges(self):
        for e in self.edgelist:
            print(e)

    def return_degree(self, dir, color):
        """Returns the degree of specific direction and color
        @param dir: direction [out/in]
        @param color: a color of an edge [act/repr/total]"""
        if dir not in ['out', 'in'] or color not in ['act', 'repr', 'total']:
            raise ValueError('Direction has to be out (or in) and color has to be act (or repr or total)!')
        if color == 'total':
            non_self_edges = getattr(self, dir + '_m')
            # non_self_edges_sum = getattr(self, dir + '_act') + getattr(self, dir + '_repr')
            # print(f'{dir} for act+repr is {non_self_edges_sum} and {dir}_m {non_self_edges}')
            self_edges = getattr(self, 'self_act') + getattr(self, 'self_repr')
        else:
            non_self_edges = getattr(self, dir + '_' + color)
            self_edges = getattr(self, 'self_' + color)
        return non_self_edges + self_edges

    def __init__(self, name):
        self.in_m = 0
        self.in_repr = 0
        self.in_act = 0
        self.out_repr = 0
        self.out_act = 0
        self.out_m = 0
        self.self_act = 0
        self.self_repr = 0
        self.edgelist = []
        self.out_act_neighbors = []
        self.out_repr_neighbors = []
        self.in_act_neighbors = []
        self.in_repr_neighbors = []
        self.name = name
        self.bconnect = {}
        self.state = None  # the state (A/R) has not be assigned
        self.R_j_estimate = None,
        self.active_out_m = 0,
        self.active_in_m = 0

    def addEdge(self, edge):
        if edge.source != self and edge.target != self:
            raise ValueError("Edge to be added have no self source nor self target!")
        self.edgelist.append(edge)
        if edge.mode != 'activation' and \
                edge.mode != 'repression' and \
                edge.mode != 'unknown':
            raise ValueError("Not supported mode of interaction!")
        if edge.source == self and edge.target != self:
            self.out_m += 1
            self.active_out_m += 1
            if edge.mode == 'activation':
                self.out_act += 1
                self.out_act_neighbors.append(edge.target)
            elif edge.mode == 'repression':
                self.out_repr += 1
                self.out_repr_neighbors.append(edge.target)
        if edge.target == self and edge.source != self:
            self.in_m += 1
            self.active_in_m += 1  # self loops are not counted in active_in_m and active_out_m
            if edge.mode == 'activation':
                self.in_act += 1
                self.in_act_neighbors.append(edge.source)
            elif edge.mode == 'repression':
                self.in_repr += 1
                self.in_repr_neighbors.append(edge.source)
        if edge.target == self and edge.source == self:
            if edge.mode == 'activation':
                self.self_act += 1
            elif edge.mode == 'repression':
                self.self_repr += 1

    def num_cons_edges(self, color, model):
        if color == 'A':
            if model in ['SC', 'SDC']:
                return self.out_act
            elif model in ['TC', 'TDC']:
                return self.in_act
        elif color == 'R':
            if model in ['SC', 'SDC']:
                return self.out_repr
            elif model in ['TC', 'TDC']:
                return self.in_repr


    def try_set_neighbors(self, params, C_node):
        #self_act_neighbor = [self] if self.self_act > 0 else []
        #self_repr_neighbor = [self] if self.self_repr > 0 else []  # self.state is necessary set
        # neighbors = self.in_act_neighbors + self.in_repr_neighbors+self.out_act_neighbors+self.out_repr_neighbors+ self_act_neighbor + self_repr_neighbor  # get N(self)
        in_neighbors = self.in_act_neighbors + self.in_repr_neighbors # N^in(self)
        out_neighbors = self.out_act_neighbors+self.out_repr_neighbors # get N^out(self)
        dirs = ['out'] * len(in_neighbors) + ['in'] * len(out_neighbors)
        N = zip(in_neighbors+out_neighbors, dirs)
        for u, dir in N:
            att = f'active_{dir}_m'
            u.__setattr__(att, getattr(u, att) - 1)  # update its neighbors
            if u.self_repr == 0 and u.self_act == 0 and u.active_in_m == 0 and u.active_out_m == 0:  # estimate independent ratio only for genes with unassigned colors, active degree 0, and no self-loops i.e. CICS
                estimate_independent_ratio(u, params, C_node)

class Edge:
    """
    A graph edge contains:
    :slot source: the source gene
    :slot target: the target gene
    :slot mode: activation/repression
    """
    __attr__ = "source", "target", "mode"

    def __str__(self):
        return ("Edge with source {0}, target {1}, mode {2}".format(self.source.name, self.target.name, self.mode))

    def __init__(self, source, target, mode):
        self.source = source
        self.target = target
        self.mode = mode


def testingGene():
    a = Gene("cex")
    b = Gene("nox")
    c = Gene("apr")
    a.addEdge(Edge(a, b, 'activation'))
    print(a)


def st_colors(dir, state_u, state_x):
    if dir == 'out':
        return f'{state_u}_{state_x}'
    elif dir == 'in':
        return f'{state_x}_{state_u}'
    else:
        raise ValueError("Only in and out dirs are allowed in st_colors")


def calculate_edge_prob(params: Dict, dir: str, e_type: str, state_u: str, state_x: str) -> float:
    param_name = 'xi_'
    if dir == 'self':
        param_name += f'{state_u}_{state_u}'
    else:
        param_name += st_colors(dir, state_u, state_x)
    if e_type == 'A':
        return params[param_name]  # prob of activation edge
    elif e_type == 'R':
        return 1 - params[param_name] # prob of repression edge
    else:
        raise ValueError("Only activation and repression edge allowed in calculate_edge_prob")



def estimate_independent_ratio(u: Gene, params: dict, C_node: List[str]) -> None:
    "Tries to estimate the ratio sigma for gene u directly. If succeeds then store R_j and the MLE color at u's attributes and updates the vertex coloring list "
    # self_repr = [u] if u.self_repr > 0 else []  # make a dummy list with u gene in case of self loops
    # self_act = [u] if u.self_act > 0 else []
    N = [u.out_repr_neighbors, u.in_repr_neighbors, u.out_act_neighbors, u.in_act_neighbors]
    edge_types = ['R', 'R', 'A', 'A']
    dirs = ['out', 'in','out', 'in']
    weights = {'A': 1,
               'R': 1}
    for neighbors, e_type, dir in zip(N, edge_types, dirs):
        for x in neighbors:
            # define the prefix for xi probability based on attr
            if x.state is None:
                raise ValueError("Trying to estimate independent ratio for non CICS! (has unassigned neighbors)")
            else:
                weights['A'] *= calculate_edge_prob(params, dir, e_type, 'A', x.state)  # append the probability mass corresponding to edges
                weights['R'] *= calculate_edge_prob(params, dir, e_type, 'R', x.state) # append the probability mass corresponding to edges
                # the order of u.state and x.state depends on the attribute type
                # for outgoing edges it is 'xi_'+u.state + '_' + x.state; for incoming: 'xi_'+x.state+'_'+u.state
                # for self loops: 'xi_'+x.state + '_' + x.state

    weights['A'] *= params['eta']# append the probability mass corresponding to vertices
    weights['R'] *= 1-params['eta']# append the probability mass corresponding to vertices
    if weights['A']/(weights['A'] + weights['R']) > 1/2:
        u.state = 'A'
        u.R_j_estimate = log10(weights['A']/(weights['A']+weights['R']))
    else:
        u.state = 'R'
        u.R_j_estimate = log10(weights['R']/(weights['A']+weights['R']))
    C_node[u.ind] = u.state # update the vertex coloring 
