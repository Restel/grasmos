import csv
import itertools
import pandas as pd
import numpy as np
import sys
from math import log10, e, log
from random import randint, choice, seed, random
from Gene import Gene, Edge, opposite_color
from lh_brute import prob_state_BNC, estimate_likelihood_brute_BNC
import pathlib
import os
from tqdm import tqdm


# logfile = f"./output/G{m}x{n}.csv" # different name from the output
# write2file(logfile, ["G_size", "H_size", "time(sec)", "num_trials_k", "timestamp"])

def write2file(logfile, text2write):
    """Use this so each line gets written to the file as soon as the
    line is ready.  Ensuring that the file is written is useful if job
    runs out of time (or space) on the cluster.  Append to the file.
    """
    f = open(logfile, "a", newline="")
    writer = csv.writer(f)
    writer.writerow(text2write)
    f.close()


def reset_R_j(genelist):
    for gene in genelist:
        gene.R_j_estimate = None
        gene.state = None
        gene.active_in_m = gene.in_m
        gene.active_out_m = gene.out_m
        # Reset bconnect as well?


def load_params(param_name):
    """Initialize the parameters for BNC model estimation for a single parameter setting.
    @param param_name: a name of a txt file with parameter values in the form:
    headers,
    xi_AA, xi_AR, xi_RA, xi_RR, eta, dataset_name
    @return params:
    """
    params = {}
    with open(param_name, 'r') as f:
        header = f.readline().rstrip().split(',')
        values = f.readline().rstrip().split(',')
        for p, val in zip(header, values):
            params[p] = float(val) if not val[-3:] == 'csv' else val
            print(p, params[p])
    return params


def createToyExamples():
    for type in range(4):
        df = toy_example(type)
        name = str(pathlib.Path().absolute()) + "/data_clean/" + "toy" + str(type) + '.csv'
        df.to_csv(name, index=False)


def createParamsMCMC():
    params_names = ['xi_A_A', 'xi_A_R', 'xi_R_A', 'xi_R_R', 'eta', 'dname']
    param_values = [0.25, 0.5, 0.75]  # 3 items, so 27 combinations
    dname_values = ['thrust_mouse.csv']

    probabilities = list(itertools.product(param_values, repeat=5))

    # Get useful lengths
    n = len(probabilities)
    k = len(dname_values)

    # The (g,) syntax converts g to a tuple, which may then be added to
    # the tuple of probabilities p
    all_params = [p + (g,) for p, g in zip(probabilities, int(n / k) * dname_values)]
    # Produce a list of dictionaries of parameters
    l = [{key: val for (key, val) in zip(params_names, ele)} for ele in all_params]
    path = pathlib.Path().absolute()
    for f in l:
        name = f"{path}/parameters/{f['xi_A_A']}_{f['xi_A_R']}_{f['xi_R_A']}_{f['xi_R_R']}_{f['eta']}_{f['dname'][:-4]}.txt"
        with open(name, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=params_names)
            writer.writeheader()
            writer.writerows([f])


def SaveParamsMCMC():

    param_value_coarse = [0.25, 0.5, 0.75]
    p = {'xi_AA': param_value_coarse,
         'xi_AR': param_value_coarse,
         'xi_RA': param_value_coarse,
         'xi_RR': param_value_coarse,
         'eta': param_value_coarse,
         'dname': ('df_regulon.csv',)}
    createParamsMCMC(p)
    p = {'xi_AA': param_value_coarse,
         'xi_AR': param_value_coarse,
         'xi_RA': param_value_coarse,
         'xi_RR': param_value_coarse,
         'eta': param_value_coarse,
         'dname': ('df_subti.csv',)}

    createParamsMCMC(p)
    p = {'xi_AA': [0.7, 0.75, 0.8],
         'xi_AR': [0.7, 0.75, 0.8],
         'xi_RA': [0.1, 0.15, 0.2, 0.25],
         'xi_RR': [0.1, 0.15, 0.2, 0.25],
         'eta': [0.4, 0.5, 0.6],
         'dname': ('df_regulon.csv',)}
    createParamsMCMC(p)
    #SUbti
    p = {'xi_AA': [0.9, 0.95, 0.99],
            'xi_AR': [0.9, 0.95, 0.99],
            'xi_RA': [0.05, 0.1, 0.15, 0.2],
            'xi_RR': [0.05, 0.1, 0.15, 0.2],
            'eta': [0.4, 0.5, 0.6],
            'dname': ('df_subti.csv',)}
    createParamsMCMC(p)



def createParamsMCMC(params):
    params_names = ['xi_A_A', 'xi_R_R', 'xi_A_R', 'xi_R_A', 'eta', 'dname']
    xi_AA_val = params['xi_AA']
    xi_AR_val = params['xi_AR']
    xi_RA_val = params['xi_RA']
    xi_RR_val = params['xi_RR']
    eta_val = params['eta']
    dname = params['dname']
    all_params = list(itertools.product(xi_AA_val,xi_RR_val,xi_AR_val,xi_RA_val, eta_val, dname))
    # Produce a list of dictionaries of parameters
    l = [{key: val for (key, val) in zip(params_names, ele)} for ele in all_params]
    path = pathlib.Path().absolute()
    for f in l:
        name = f"{path}/parameters/{f['xi_A_A']}_{f['xi_A_R']}_{f['xi_R_A']}_{f['xi_R_R']}_{f['eta']}_{f['dname'][:-4]}.txt"
        with open(name, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=params_names)
            writer.writeheader()
            writer.writerows([f])

def CreateExamplesGradDescent():
    true_params = {'xi_A_A': 0.3,
                   'xi_A_R': 0.2,
                   'xi_R_A': 0.7,
                   'xi_R_R': 0.9,
                   'eta': 0.6,
                   'lamda': 0.05,
                   'final_precision': 0.01
                   }
    for ind in range(5):  # make 5 datasets of the same parameter set
        # create columns for pd dataframe
        source = []
        target = []
        mode = []
        # assign roles to genes
        colors = []
        for _ in range(200):
            if random() < true_params['eta']:
                colors.append('A')
            else:
                colors.append('R')
        # create 40 genetic interactions based on gene's color assignment
        for _ in range(20):
            ind_s = randint(0, 199)
            ind_t = randint(0, 199)  # create an edge
            source.append("v" + str(ind_s))
            target.append("v" + str(ind_t))
            pair = f'{colors[ind_s]}_{colors[ind_t]}'  # color the edge
            if random() < true_params['xi_' + pair]:
                mode.append("activation")
            else:
                mode.append("repression")
        df = pd.DataFrame(dict(source=source, target=target, mode=mode))
        name = str(pathlib.Path().absolute()) + "/data_clean/" + "examplesGD" + str(ind) + '.csv'
        df.to_csv(name, index=False)


def toy_example(type):
    seed(123)
    source = []
    target = []
    mode = []
    if type == 0:
        source = ['v1', 'v2', 'v2', 'v2']
        target = ['v1', 'v1', 'v1', 'v1']
        mode = ['repression', 'activation', 'activation', 'repression']
    if type == 1:  # all links are activators
        for i in range(6):
            for j in range(6):
                source.append('v' + str(i))
                mode.append('activation')
                ind = randint(0, 5)
                target.append('v' + str(ind))
    if type == 2:  # all links are repression
        for i in range(6):
            for j in range(6):
                source.append('v' + str(i))
                mode.append('repression')
                ind = randint(0, 5)
                target.append('v' + str(ind))
    if type == 3:  # random activation or repression
        for i in range(6):
            for j in range(6):
                source.append('v' + str(i))
                mode_ch = choice(['activation', 'repression'])
                ind = randint(0, 5)
                target.append('v' + str(ind))
                mode.append(mode_ch)
    df = pd.DataFrame(dict(source=source, target=target, mode=mode))
    return df


def initialize_toy_graph(toy_examples, k):
    """Function to initialize graphs for toy examples of MCMC sampling. toy_examples is a list of strings, where each string is in form 'source_i,target_i,mode_i;source_i+1,target_i+1,mode_i+1; ...' for i in [1, m]
    @param toy_examples: list of string toy examples
    @param k: id of toy example
    """
    ex = toy_examples[k]
    strings = ex.split(';')
    toy_df = []
    for edge in strings:
        toy_df.append(edge.split(','))
    df = pd.DataFrame(toy_df, columns=['source', 'target', 'mode'])
    genelist, genedict = initialize_graph(df)
    return genelist, genedict, df


def initialize_params(toy_params, k):
    """Function to initialize parameter settings for toy examples of MCMC sampling. toy_params is a list of strings, where the first row is the string of parameter names, followed by strings of parameter values
        @param toy_params: list of parameter values for each toy setting
        @param k: id of toy example
        """
    params_names = toy_params[0].split(';')  # get parameter names
    params_values = toy_params[k].split(';')  #
    params = {key: float(val) for (key, val) in zip(params_names, params_values)}
    return params


def load_data(dataname):
    pathname = str(pathlib.Path().absolute()) + "/data_clean/"
    data = pd.read_csv(pathname + dataname + ".csv")
    return data


class Dataframe:
    def __init__(self, dataname):
        self.dataname = dataname
        if isinstance(dataname, str):
            self.df = load_data(dataname)  # constructor loading a datafile by its name
        elif isinstance(dataname, pd.DataFrame):
            self.df = dataname  # constructor loading a datafile from pandas dataframe
        self.E_A = self.df.query("mode == 'activation'").shape[0]  # number of activation edges
        self.E_R = self.df.query("mode == 'repression'").shape[0]  # number of repression edges


def ER_likelihood(E_A, E_R):
    """Negative loglikehood for Erdos-Renyi model"""
    prob_a = min(E_A / (E_A + E_R), 0.99999999) # not let prob_a be 1
    return -1 * (E_A * log10(prob_a) + E_R * log10(1 - prob_a))


def exponent_repr(node, model):  # exponent for inconsistent edges when the gene is activator
    if model == 'SC' or model == 'SDC':
        exp = node.return_degree('out', 'repr')  # node.out_repr + node.self_R
    elif model == 'TC' or model == 'TDC':
        exp = node.return_degree('in', 'repr')  # node.in_repr + node.self_R
    else:
        raise ValueError("model {0} is not supported".format(model))
    return exp


def exponent_act(node, model):
    if model == 'SC' or model == 'SDC':
        exp = node.return_degree('out',
                                 'act')  # in_act include outgoing edges excluding self loops node.out_act + node.self_A
    elif model == 'TC' or model == 'TDC':
        exp = node.return_degree('in', 'act')  # in_act include incoming edges excluding self loops
    else:
        raise ValueError("model {0} is not supported".format(model))
    return exp


def edges_cons(node, color, model):
    if color == 'A':
        return (exponent_act(node, model))
    if color == 'R':
        return (exponent_repr(node, model))


def edges_incons(node, color, model):
    if color == 'A':
        return (exponent_repr(node, model))
    if color == 'R':
        return (exponent_act(node, model))


def logistic_fun(node, model, gamma):
    if model == 'SDC':
        return 1 - 1 / 2 * e ** (-gamma * node.return_degree('out',
                                                             'total'))  # we added self connections to correspond to the test suit.
    elif model == 'TDC':
        return 1 - 1 / 2 * e ** (-gamma * node.return_degree('in',
                                                             'total'))  # we added self connections to correspond to the test suit. node.in_m + self_connections


def xi_act(xi_A, node, model, map_A=logistic_fun, **kwargs):
    if model == 'SC' or model == 'TC':
        prob = xi_A
    elif model == 'SDC' or 'TDC':
        prob = map_A(node, model, **kwargs)
    else:
        raise ValueError("model {0} is not supported".format(model))
    return abs(prob - 10 ** -10)


def xi_repr(xi_R, node, model, map_R=logistic_fun, **kwargs):
    if model == 'SC' or model == 'TC':
        prob = xi_R
    elif model == 'SDC' or 'TDC':
        prob = map_R(node, model, **kwargs)
    else:
        raise ValueError("model {0} is not supported".format(model))
    return abs(prob - 10 ** -10)


def xi(xi_dict, gene, color, model, **kwargs):
    if color == "A":
        return xi_act(xi_dict['A'], gene, model, **kwargs)
    elif color == "R":
        return xi_act(xi_dict['R'], gene, model, **kwargs)


def sigma(node, xi_A, xi_R, model, eta, **kwargs):
    """Estimates the maximum ratio of two partition function: Z_j and Z_j-1 and dynamically chooses the next color in
    the chain of products """
    cons_A = xi_act(xi_A, node, model, **kwargs)  # estimate the consistency probability between the activation
    # node and the activation edge for the specific node under the specific model
    cons_R = xi_repr(xi_R, node, model, **kwargs)  # estimate the consistency probability for the node under the model

    try:
        alpha_a = (cons_R / (1 - cons_A)) ** exponent_repr(node, model) * ((1 - cons_R) / cons_A) ** exponent_act(node,
                                                                                                                  model) * (
                          1 - eta) / eta
    except OverflowError:
        return 0, 'R'
    if alpha_a == 0:
        return 0, 'A'  # -log10(1) = 0
    elif alpha_a < 1 / alpha_a:
        alpha = alpha_a
        color = 'A'
    else:
        alpha = 1 / alpha_a
        color = 'R'
    return - log10(1 + alpha), color


def Z_all_act(genelist, model, xi_A, eta, **kwargs):
    """Function to compute the conditional probability of a set of colored edges, given the graph topology, a coloring C
    (with all genes set as activators), xi and eta, i.e. returns Pr(E_C|G,C, xi, eta) * Pr(C) """
    log_lh = 0
    for gene in genelist:
        prob_act = xi_act(xi_A, gene, model, **kwargs)
        exp_act = exponent_act(gene, model)
        exp_repr = exponent_repr(gene, model)
        log_lh = log_lh + exp_act * log10(prob_act) + exp_repr * log10(1 - prob_act)
    return log_lh + len(genelist) * log10(eta)


def Z_all(genelist, model, xi_A, xi_R, eta, coloring, **kwargs):
    """Function to compute the conditional probability of a set of colored edges, given the graph topology, an arbitrary node coloring C, consistency probabilities, and eta, i.e. returns Pr(E_C|G,C, xi, eta) * Pr(C) """
    log_lh = 0
    xi_dict = {'A': xi_A,  # model parameters
               'R': xi_R}
    q_dict = {'A': eta,  # model parameters
              'R': 1 - eta}
    for i in range(len(genelist)):
        gene = genelist[i]
        color = coloring[i]
        prob_cons = xi(xi_dict, gene, color, model, **kwargs)
        exp_cons = edges_cons(gene, color, model)
        exp_incons = edges_incons(gene, color, model)
        log_lh = log_lh + exp_cons * log10(prob_cons) + exp_incons * log10(1 - prob_cons) + log10(q_dict[color])
    return log_lh


def computeZ(genelist, model, xi_A, xi_R, eta, **kwargs):
    """Function to compute a partition function for defined consistency probability xi, and activation probability
    eta """
    log_product = 0
    coloring_nodes = ''
    for j in range(0, len(genelist)):
        log_add, color = sigma(genelist[j], xi_A, xi_R, model, eta, **kwargs)
        # print(f'log for {j}th vertex sigma is {log_add}')
        log_product += log_add
        coloring_nodes += color
    Z_all_log = Z_all(genelist, model, xi_A, xi_R, eta, coloring_nodes,
                      **kwargs)
    # print(f'Z for all nodes colored: {Z_all_log}')
    # print(f'MLE node coloring:{coloring_nodes}')
    print(Z_all_log)
    return Z_all_log - log_product, coloring_nodes


def calculate_transition_prob(params, gene, C_node):
    """Calculate the transition probability from a state where gene is activator, to a state where it is a repressor given the colors of other nodes as specified in C_node.
    params is a dict with the following keys: 'xi_A_A,  'xi_R_A', 'xi_A_R', 'xi_R_R' and values in [0,1] """
    gene_neighbors = gene.bconnect  # O(1)
    to_activators = gene_neighbors['out_A_A'] * log10(params['xi_R_A'] / params['xi_A_A']) + gene_neighbors[
        'out_R_A'] * log10((1 - params['xi_R_A']) / (1 - params[
        'xi_A_A']))  # problog corresponding to the ratio change related to the outgoing edges from the gene to activators
    to_repressors = gene_neighbors['out_A_R'] * log10(params['xi_R_R'] / params['xi_A_R']) + gene_neighbors[
        'out_R_R'] * log10((1 - params['xi_R_R']) / (1 - params['xi_A_R']))
    from_activators = gene_neighbors['in_A_A'] * log10(params['xi_A_R'] / params['xi_A_A']) + gene_neighbors[
        'in_R_A'] * log10((1 - params['xi_A_R']) / (1 - params['xi_A_A']))
    from_repressors = gene_neighbors['in_A_R'] * log10(params['xi_R_R'] / params['xi_R_A']) + gene_neighbors[
        'in_R_R'] * log10((1 - params['xi_R_R']) / (1 - params['xi_R_A']))
    vertex_color = log10((1 - params['eta']) / params['eta'])
    self_activation = gene_neighbors['self_A'] * log10((params['xi_R_R'] / params['xi_A_A']))
    self_repression = gene_neighbors['self_R'] * log10((1 - params['xi_R_R']) / (1 - params['xi_A_A']))
    log_prob = to_activators + to_repressors + from_activators + from_repressors + vertex_color + self_activation + self_repression
    if log_prob > 308:  # to avoid OverflowError
        output = sys.float_info.max
    elif log_prob < -300:  # to avoid ZeroDivision Error
        output = sys.float_info.min
    else:
        output = 10 ** log_prob
    return output


def MC_sample(C_node, j, genelist, T, params):
    """Function for sampling a node coloring according to a stationary distribution defined by edge coloring and BNC model
    @param C_node: a partial node_coloring with first j vertices colored, and n-j remaining vertices uncolored, e.g. ['A', 'A', None, None, etc]
    @param j: number of colored nodes
    @param genelist: a list of Gene objects
    @param T: mixing time, the number of steps to run MC
    @param params: BNC parameters
    @return C_node_new: return the resulting node_coloring after T steps
    """
    pi = C_node.copy()
    sample_indices = []
    for ind, elem in enumerate(pi):
        if elem is None:
            pi[ind] = choice(['A', 'R'])
            sample_indices.append(ind)
    for gene in genelist:
        gene.estimate_bi_node_connections(pi)  # estimate the initial bconnections for all genes
    for t in range(T):
        # if t % 10000 == 0:
        #     print(f'Done for {t} iteration, {j} vertex')
        # l = np.random.choice(np.arange(j, n))  # choose an uncolored remaining gene u.a.r
        # np.arrange(j,n) = array([j, ...,n-1]) and random.randint(j, n) includes both j and n as valid choices
        # l = randint(j, n - 1)  # choose an uncolored remaining gene u.a.r
        # list of gene indices + boolean mask -> give the array of the unprocessed lists for random.choice()
        # if its empty return
        l = choice(sample_indices)
        gene = genelist[l]
        prob_A_R = calculate_transition_prob(params, gene, pi)
        random_number = random()
        prob = min(1, prob_A_R) if pi[l] == 'A' else min(1, 1 / prob_A_R)
        if random_number < prob:
            pi[l] = opposite_color(pi[l])
            gene.update_connections(pi[l])
    return pi


def Z_all_MC(C_node, df, genedict, params):
    """Function for computing the likelihood of an edge coloring given the node coloring under BNC model. Returns a log likelihood
    @param C_node: a full node coloring
    @param df: a dataframe where each row is an edge with str source target type
    @param genelist: a list of Gene objects
    @param gene_dict: a dict of gene names and their corresponding indices
    @param params: dictionary with model parameters
    """
    log_prod = 0
    for ind in df.index:
        target_ind = genedict[df['target'][ind]]
        target_color = C_node[target_ind]
        source_ind = genedict[df['source'][ind]]
        source_color = C_node[source_ind]
        if df['mode'][ind] == 'activation':
            log_prod += log10(params['xi_' + source_color + '_' + target_color])
        elif df['mode'][ind] == 'repression':
            log_prod += log10(1 - params['xi_' + source_color + '_' + target_color])
        else:
            raise ValueError(f"Unsuported interaction mode {df['mode'][ind]}")
    n_act = C_node.count('A')
    n_repr = len(C_node) - n_act
    log_prod += n_act * log10(params['eta']) + n_repr * log10(1 - params['eta'])
    # print(f'The likelihood for MLE color partition {C_node}: {10 ** log_prod}')
    return log_prod


def derivative(p, edge_df, C_node, params, genedict):
    """Estimates the partial derivative of the log likelihood by the parameter"""
    # df = df[df['mode'] != 'unknown'] we are assuming that all unknown edges have been filtered out
    if p == 'eta':
        n_act = C_node.count('A')
        n_repr = len(C_node) - n_act
        d = n_act / (params[p] * log(10)) - n_repr / (log(10) - params[p] * log(10))
    else:  # derivative by one of the edge probabilities
        e_act = 0
        e_repr = 0
        param_source, param_target = p[3], p[5]
        for ind, row in edge_df.iterrows():
            source_ind = genedict[row['source']]
            target_ind = genedict[row['target']]
            source_color = C_node[source_ind]
            target_color = C_node[target_ind]
            if source_color == param_source and target_color == param_target:
                if row['mode'] == 'activation':
                    e_act += 1
                elif row['mode'] == 'repression':
                    e_repr += 1
                else:
                    raise ValueError("Unsupported edge type")
        d = e_act / (params[p] * log(10)) - e_repr / (log(10) - params[p] * log(10))
    return d


def estimate_gradient(C_node, params, edge_df, genedict):
    gradient = []
    p_names = ['xi_A_A', 'xi_A_R', 'xi_R_A', 'xi_R_R', 'eta']
    for p in p_names:
        gradient.append(derivative(p, edge_df, C_node, params, genedict))
    return gradient


def update_params(params, gradient):
    """Function to updates the current values of activation probabilities in params dict"""
    p_names = ['xi_A_A', 'xi_A_R', 'xi_R_A', 'xi_R_R', 'eta']
    for p_name, p_val in zip(p_names, gradient):
        params[p_name] = min(0.999, params[p_name] + p_val * params['lambda'])


def check_convergence(cur_params, prev_params, eps):
    """Function to updates the current values of activation probabilities in params dict"""
    p_names = ['xi_A_A', 'xi_A_R', 'xi_R_A', 'xi_R_R', 'eta']
    delta = []
    for p_name in p_names:
        delta.append(abs(cur_params[p_name] - prev_params[p_name]))
    cur_precision = max(delta)
    if cur_precision < eps:
        converged = True
    else:
        converged = False
    return converged, cur_precision


def computeZBinode(genelist, df, genedict, params, T, S, eps):
    """A function to compute a partition function for Binode consistent model given the parameters defined in params and the edge coloring
    @param genelist: a list of Gene objects ordered according to their indices
    @param genedict: a dict with genes' names and their indices
    @param df: a dataframe of edges
    @param params: a dictionary with 4 activation edge probabilities, node activation probability, gamma (if applicable)
    @param eps: convergence error for testing MC convergence
    @param T: number of time steps for each MC
    @param S: number of samples for MC
    """
    # if 'dataset' in params.keys():
    #     filename = f"/results/MC_convergence_test_{params['dataset']}.txt"
    # else:
    #     filename = '/results/MC_convergence_test.txt'
    C_node = [None] * len(genelist)  # a vector of assigned colors to genes -> an MLE coloring
    log_sum_R = 0
    # sort the genes by their total degree (in+out, both act and repr)
    # update the dictionary with the current indices
    gradient = [0] * 5
    for j in tqdm(range(0, len(genelist))):
        # print(f'Done for {j} vertex')
        # create a list of unprocessed genes, if empty then continue (proceed to the next gene calculation_
        gene_j = genelist[j]
        if gene_j.state is not None:  # the gene j and its corresponding sigma_j ratio has already been processed as a part of the independent island calculation
            log_sum_R += gene_j.R_j_estimate
            C_node[j] = gene_j.state  # make sure that j is the right index of the gene_j
            print(f'{j} vertex CICS, Rj {10 ** gene_j.R_j_estimate}')  # will add a newline for the next vertex
            continue
        sum = 0
        R_arr = [None] * 10  # array with Zj/Z_[j-1] estimates
        for s in tqdm(range(S)):
            C_sampled = MC_sample(C_node, j, genelist, T, params)
            if C_sampled[j] == 'A':
                sum += 1
            R_arr[s % 10] = sum / (s + 1)  # collecting last 10 R estimates to check the convergence
            # print(f'Vertex {j}, iteration {s}')
            # sample_gradient = estimate_gradient(C_node, params, df,
            #                                     genedict)  # estimate the log gradient by the sampled node coloring
            print(sum / (s + 1), end=',')  # print current estimate of R_j for diagnostics
            # if j == 0:
            #     gradient[:] = [x + y for x, y in zip(gradient, sample_gradient)]  # add the gradient of the sample
            if s > 10 and all(abs(np.diff(R_arr)) < eps):
                break
        # gradient = gradient/(s+1)  # average by the actual number of samples taken (<= S)
        R_j = sum / (s + 1)
        if R_j > 1 / 2:
            C_node[j] = 'A'
            gene_j.state = 'A'
            # gene_j.R_j_estimate = log10(R_j)
            log_sum_R += log10(R_j)
        else:
            C_node[j] = 'R'
            gene_j.state = 'R'
            # gene_j.R_j_estimate = log10(1-R_j)
            log_sum_R += log10(1 - R_j)
        print(j, R_j)  # will add a newline for the next vertex
        # if j == 0:
        #     gradient[:] = [x / (s+1) for x in gradient]  # average by the total number of samples taken (<= S)
        gene_j.try_set_neighbors(params,
                                 C_node)  # try setting up the colors and the ratios as a part of independent island calculation # walk over each connection v_u of v_j, and check its neighbors, if they are only connected to colored genes -> estimate sigma_j in a closed form from parameters 
        # try_set_neighbors also updates input argument C_node node coloring
    Zn = Z_all_MC(C_node, df, genedict, params)
    # write Z_all for MLE_color
    print(Zn)
    # write the MLE_color
    MLE_node = ''.join(C_node)
    print(MLE_node)
    return Zn - log_sum_R, gradient, MLE_node  # save the result into a new file + the log portion from all the independent subgraphs (log_independent_subraphs)


def compute_likelihood(genelist, model, **kwargs):
    """Negative loglikehood for SC, TC, SDC, TDC models. For SDC and TDC gamma parameter is needed"""
    gradient = [0, 0, 0, 0, 0]
    MLE_node_coloring = None
    Z = None
    if model == 'BNC':
        Z, gradient, MLE_node_coloring = computeZBinode(genelist, **kwargs)
    elif model in ['SC', 'TC', 'SDC', 'TDC']:
        Z, MLE_node_coloring = computeZ(genelist, model, **kwargs)
    return -1 * Z, gradient, MLE_node_coloring


def estimate_likelihood(model, df, gene_list, **kwargs):
    lh = None
    gradient = None
    node_coloring = None
    if model == 'ER':
        lh = ER_likelihood(df.E_A, df.E_R)
    else:
        lh, gradient, node_coloring = compute_likelihood(gene_list, model, **kwargs)
    return lh, gradient, node_coloring  # FIX THIS to full precision?


def gene_names(df):
    """Returns a list of unique gene names from a dataframe with edges in the form [source, target, mode]"""
    names = pd.concat([df['source'], df['target']]).unique()
    return names


def initialize_genes(names, genelist, genedict):
    i = 0
    for name in names:
        gene = Gene(str(name))
        gene.ind = i
        gene.active_in_m = gene.in_m
        gene.active_out_m = gene.out_m
        genedict[name] = i
        genelist.append(gene)
        i += 1


def initalize_edges(df, genelist, genedict):
    for ind in df.index:
        target_ind = genedict[df['target'][ind]]
        source_ind = genedict[df['source'][ind]]
        target = genelist[target_ind]
        source = genelist[source_ind]
        if source != target:  # to avoid double counting in self-loops
            source.addEdge(Edge(source, target, df['mode'][ind]))
            target.addEdge(Edge(source, target, df['mode'][ind]))
        else:
            source.addEdge(Edge(source, target, df['mode'][ind]))


def initialize_graph(df):
    names = gene_names(df)
    genelist = []
    genedict = {}
    initialize_genes(names, genelist, genedict)
    initalize_edges(df, genelist, genedict)
    return genelist, genedict


def try_append_results(csv_file, csv_columns, dict_data):
    try:
        with open(csv_file, 'a') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writerow(dict_data)
    except IOError:
        print("I/O error")


def write_headers(csv_file, headers):
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()  # write column headers


def return_param_set(model):
    """Returns a list of dictionaries - all combinations of parameters to try for a given model. SC and TC have
    parameters xi_A, xi_R, eta; SDC and TDC additionally have gamma """
    params_names = ['xi_A', 'xi_R', 'eta']
    param_values = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    if model in ['SDC', 'TDC']:
        temp = list(itertools.product(param_values, repeat=4))
        params = [{key: val for (key, val) in zip(params_names + ['gamma'], ele)} for ele in temp]
    elif model in ['SC', 'TC']:
        temp = list(itertools.product(param_values, repeat=3))
        params = [{key: val for (key, val) in zip(params_names, ele)} for ele in temp]
    else:
        raise ValueError("Model is not supported")
    return params


def estimate_stat_distro(df, params, gene_dict):
    Omega = list(itertools.product(['A', 'R'], repeat=len(gene_dict)))
    stat_distro = {}
    Z = 0
    for w in Omega:
        state_w = prob_state_BNC(w, gene_dict, df, params)
        stat_distro[w] = state_w
        Z += state_w
    for w in Omega:
        stat_distro[w] = stat_distro[w] / Z
    return stat_distro


def estimate_TV(C_set, n, stat_distro):
    Omega = list(itertools.product(['A', 'R'], repeat=n))
    TV = 0
    for w in Omega:
        # w_count = prob_state_BNC(list(w), gene_dict, df, params)
        w_count = stat_distro[w]  # implement a stationary distribution
        if w_count > 1:
            raise ValueError("Probability can not exceed 1!")
        w_obs_count = C_set[w] / C_set['s'] if w in C_set.keys() else 0
        TV += abs(w_count - w_obs_count)
    frac = (len(C_set) - 1) / 2 ** n
    return TV * 0.5, frac


def test_MC_sample(S, T, i):
    name = str(pathlib.Path().absolute()) + "/data_clean/" + f'toy{i}.csv'
    df = pd.read_csv(name)
    genelist, genedict = initialize_graph(df)
    n = len(genelist)
    d = {'xi_A_A': 0.7,
         'xi_A_R': 0.3,
         'xi_R_A': 0.6,
         'xi_R_R': 0.5,
         'eta': 0.7,
         'dataset': f'toy{i}'}
    C_node = []
    C_set = {}
    C_set['s'] = 0
    stat_distro = estimate_stat_distro(df, d, genedict)
    for s in range(S):
        C_set['s'] += 1
        C_node = MC_sample(C_node, 0, genelist, T, d)
        if tuple(C_node) in C_set.keys():
            C_set[tuple(C_node)] += 1
        else:
            C_set[tuple(C_node)] = 1
        if s % 20 == 0:
            TV, frac = estimate_TV(C_set, n, stat_distro)
            print(f'Sample {s}, TV {TV}, Percentage of states {frac}')


def main():
    out_path = str(pathlib.Path().absolute())
    datasets = ['thrust_mouse', 'df_regulon']
    headers = ['model', 'dataset', 'xi_A', 'xi_R', 'eta', 'gamma', 'lh', 'MLE_Cn']
    models = ['SDC', 'TDC', 'SC', 'TC']
    for dataname in datasets:
        csv_file = out_path + '/results/' + dataname + '.csv'
        write_headers(csv_file, headers)  # write headers per each df
        df = Dataframe(dataname=dataname)  # initialize a df
        genelist, genedict = initialize_graph(df.df)
        genelist.sort(key=lambda g: g.out_m + g.in_m, reverse=True)
        for ind, g in enumerate(genelist):
            g.ind = ind
            genedict[g.name] = ind  # updating the indices
        for mode in models:
            params = return_param_set2()
            for p in params:
                get_lh_results(dataname, mode, df, genelist, headers, csv_file, p)


def find_dual_edges(dataname):
    df = Dataframe(dataname=dataname).df
    df['id'] = df['source'] + df['target']
    duplicated = pd.concat(g for _, g in df.groupby("id") if len(g) > 1)
    duplicated['idc'] = duplicated['id'] + duplicated['mode']
    duplicated_by_color = pd.concat(g for _, g in duplicated.groupby('idc') if len(g) > 1)
    print(f'DF {dataname}, duplicated:{len(duplicated)}, percentage: {len(duplicated)/2/len(df)}, abs_duplicates:{len(duplicated_by_color)}')
    # count = {}
    # count_dual = []
    # for index,row in df.iterrows():
    #     s,t,m = row['source'], row['target'], row['mode']
    #     if (s,t) in count.keys():
    #         count[(s,t)] += 1
    #     else:
    #         count[(s, t)] = 1
    # for comb in count.keys():
    #     if count[comb] >1:
    #         count_dual.append(comb)

def get_lh_results(dataname, model, df, genelist, headers, csv_file, params):
    lh, gradient, MLE_cn = estimate_likelihood(model,
                                               df,
                                               genelist,
                                               xi_A=params['xi_A'],
                                               xi_R=params['xi_R'],
                                               eta=params['eta'],
                                               gamma=params['gamma'])
    new_row = {'model': model,
               'dataset': dataname,
               'xi_A': params['xi_A'],
               'xi_R': params['xi_R'],
               'eta': params['eta'],
               'gamma': params['gamma'],
               'lh': lh,
               'MLE_Cn': MLE_cn}
    try_append_results(csv_file, headers, new_row)
    print("Model: {0}, eta:{1}, xi_A:{2}, xi_R:{3}, gamma: {4}, lh:{5}".format(model,
                                                                               params['eta'],
                                                                               params['xi_A'],
                                                                               params['xi_R'],
                                                                               params['gamma'],
                                                                               lh))


def create_gene_order_map(dname):
    name = str(pathlib.Path().absolute()) + "/data_clean/" + dname
    df = pd.read_csv(name)
    df = df[df['mode'] != 'unknown']  # ignore the unknown edges
    df = df.reset_index()
    genelist, genedict = initialize_graph(df)
    genelist.sort(key=lambda g: g.out_m + g.in_m, reverse=True)
    for ind, g in enumerate(genelist):
        g.ind = ind
        genedict[g.name] = ind  # updating the indices
    w = csv.writer(open(f"gene_map_{dname}", "w"))
    for key, val in genedict.items():
        # write every key and value to file
        w.writerow([key, val])


def get_undefined_genes(file_name):
    # filenames = '0.75_0.75_0.25_0.25_0.5_df_regulon.txt' or '0.75_0.75_0.25_0.25_0.75_thrust_mouse.txt'
    #file_name = 'subti_refined/0.99_0.99_0.2_0.15_0.5_df_subti.txt'
    #file_name = 'regulon_refined_duplicates/0.7_0.8_0.25_0.1_0.5_df_regulon.txt'
    name = str(pathlib.Path().absolute()) + "/output/" + file_name
    with open(name) as f:
        lines = f.readlines()
    lines = lines[6:]  # remove input probs
    RA_string = lines[-4].rstrip()
    lines = lines[:len(lines) - 5]  # remove other statistis, len(lines) == 1922 for regulon and 1858 for mouse
    probabilities = []
    for l in lines:
        elems = l.rstrip().split(' ')
        ele = float(elems[-1])
        probabilities.append(ele)
    #print(probabilities)
    U_string = ['U' if ele == 0.5 else None for ele in probabilities]
    final_string = [u_state if u_state else regulator_state for regulator_state, u_state in zip(RA_string, U_string)]
    return ''.join(final_string)



def initialize_graph_instance(dname, syn=False):
    if syn:
        name = f"{pathlib.Path().absolute()}/{dname}"
    else:
        name = f"{pathlib.Path().absolute()}/data_clean/{dname}.csv"
    df = pd.read_csv(name)  # sys.argv[0] is the name of the python script
    df = df[df['mode'] != 'unknown']  # ignore the unknown edges
    df = df.reset_index()  # make sure indexes pair with number of rows
    genelist, genedict = initialize_graph(df)
    genelist.sort(key=lambda g: g.out_m + g.in_m, reverse=True)
    for ind, g in enumerate(genelist):
        g.ind = ind
        genedict[g.name] = ind  # updating the indices
    return df, genelist, genedict