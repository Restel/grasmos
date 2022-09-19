import pandas as pd
import numpy as np
from math import log10, e
from Gene import Gene, Edge
import pathlib
import itertools


def estimate_likelihood_brute(model, gene_list, **kwargs):
    Omega = list(itertools.product([0, 1], repeat=len(gene_list)))
    Z = 0
    for state in Omega:
        Z += prob_state(state, model, gene_list, **kwargs)
    return -log10(Z)


def prob_state(state, model, gene_list, eta, xi_A, xi_R, gamma = None):
    """A function to estimate the probability of a node coloring state given colored edges. Works only for one-end-point models such as SC, TC, SDC, TDC
    """
    prob_omega = 1
    xi_A = xi_A
    xi_R = xi_R
    for ind in range(0, len(gene_list)):
        gene = gene_list[ind]
        if state[ind] == 1:
            if model == 'SC':
                prob_gene = xi_A ** gene.return_degree('out', 'act') * (1-xi_A) ** gene.return_degree('out', 'repr')
            elif model == 'TC':
                prob_gene = xi_A ** gene.return_degree('in', 'act') * (1-xi_A) ** gene.return_degree('in', 'repr')
            elif model == 'SDC':
                prob_xi = 1 - 1/2 * e**(-gamma*gene.return_degree('out', 'total'))
                prob_gene = prob_xi ** gene.return_degree('out', 'act') * (1 - prob_xi) ** gene.return_degree('out', 'repr')
            elif model == 'TDC':
                prob_xi = 1 - 1/2 * e ** (-gamma*gene.return_degree('in', 'total'))
                prob_gene = prob_xi ** gene.return_degree('in', 'act') * (1 - prob_xi) ** gene.return_degree('in', 'repr')
        else:
            if model == 'SC':
                prob_gene = xi_R ** gene.return_degree('out', 'repr') * (1-xi_R) ** gene.return_degree('out', 'act')
            elif model == 'TC':
                prob_gene = xi_R ** gene.return_degree('in', 'repr') * (1-xi_R) ** gene.return_degree('in', 'act')
            elif model == 'SDC':
                prob_xi = 1 - 1/2 * e**(-gamma*gene.return_degree('out', 'total'))
                prob_gene = prob_xi ** gene.return_degree('out', 'repr') * (1-prob_xi) ** gene.return_degree('out', 'act')
            elif model == 'TDC':
                prob_xi = 1 - 1/2 * e ** (-gamma*gene.return_degree('in', 'total'))
                prob_gene = prob_xi ** gene.return_degree('in', 'repr') * (1-prob_xi)**gene.return_degree('in', 'act')
        prob_omega *= prob_gene
    n_act = sum(state)
    n_repr = len(gene_list) - n_act
    return prob_omega * eta ** n_act * (1-eta)**n_repr


def estimate_likelihood_brute_BNC(gene_dict, df, d):
    Omega = list(itertools.product(['A', 'R'], repeat=len(gene_dict)))
    Z = 0
    for state in Omega:
        prob_state = prob_state_BNC(state, gene_dict, df, d)
        #print(f'Probability of state {state} is {prob_state}')
        Z += prob_state
    return -log10(Z)


def prob_state_BNC(state, gene_dict, df, d):
    prob_omega = 1
    df = df[df['mode'] != 'unknown']
    for ind in df.index:
        source_node = df['source'][ind]
        source_ind = gene_dict[source_node]
        target_node = df['target'][ind]
        target_ind = gene_dict[target_node]
        source_color = state[source_ind]
        target_color = state[target_ind]
        prob_activation_edge = d['xi_'+ source_color + '_' + target_color]
        if df['mode'][ind] == 'activation':
            prob_edge = prob_activation_edge
        else:
            prob_edge = 1 - prob_activation_edge
        prob_omega *= prob_edge
    n_act = state.count('A')
    n_repr = len(state) - n_act
    return prob_omega * d['eta'] ** n_act * (1-d['eta'])**n_repr

