from scripts.reaction_class import Reaction, ReactionSystem
from scripts.grid_class import GridParms
from scripts.tree_class import Tree
from scripts.initial_condition_class import InitialCondition
import sympy as sp
import numpy as np
import os

class Model:

    def __init__(self, _species):

        self.reactions = []
        self.species = _species

    """
    Different function to add single or multiple reactions to our model
    """
    def add_reaction(self, reactants, products, propensities):
        
        num_symbols = len(self.species)
        reactants = reactants * (-1)
        eq_sp = reactants + products

        nu_vec = np.zeros(num_symbols)
        for i, sym in enumerate(self.species):
            nu_vec[i] = eq_sp.coeff(sym)

        prop_dict = {}
        for key, value in list(propensities.items()):
            for i, sym in enumerate(self.species):
                if(key == sym):
                    prop_dict[i] = sp.lambdify(sym, value)
        
        self.reactions.append(Reaction(prop_dict, nu_vec))

    def add_reactions(self, reactants_list, products_list, propensities_list):
        for (reactants, products, propensities) in zip(reactants_list, products_list, propensities_list):
            self.add_reaction(reactants, products, propensities)
    
    def generate_reaction_system(self):
        self.reaction_system = ReactionSystem(self.reactions, self.species)



class Partitioning:     # maybe better to call it tree?

    def __init__(self, _partition, _r, _model):
        
        self.r = _r
        self.model = _model
        self.partition = _partition
        for i, sym in enumerate(self.model.species):
            self.partition = self.partition.replace(str(sym), str(i))

    def add_grid_params(self, n, binsize, liml):

        self.grid = GridParms(n, binsize, liml)

    def generate_tree(self):

        self.tree = Tree(self.partition, self.grid)
        self.tree.initialize(self.model.reaction_system, self.r)

    def generate_initial_condition(self, n_basisfunctions):

        self.initial_conditions = InitialCondition(self.tree, n_basisfunctions)

    

def run(partitioning, input, output, snapshot, tau, tfinal, substeps, method):
    partitioning.tree.write()
    cmd = f'bin/hierarchical-cme -i {input} -o {output} -s {snapshot} -t {tau} -f {tfinal} -n {substeps} -m {method}'
    os.system(cmd)