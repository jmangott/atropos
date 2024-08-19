from scripts.reaction_class import Reaction, ReactionSystem
from scripts.grid_class import GridParms
from scripts.tree_class import Tree
from scripts.initial_condition_class import InitialCondition
from scripts.index_functions import incrVecIndex
import sympy as sp
import numpy as np
import os
import sys

class Model:

    def __init__(self, _species):

        self.reactions = []
        self.species = _species

    """
    Different function to add single or multiple reactions to our model
    """
    def add_reaction(self, reactants, products, propensities):

        num_symbols = len(self.species)
        eq_sp = products - reactants

        nu_vec = np.zeros(num_symbols)
        for i, sym in enumerate(self.species):
            nu_vec[i] = eq_sp.coeff(sym)

        prop_dict = {}
        # Test if we only have coefficient as variable, if so, generate propensity in non factorised form
        if type(propensities) == int or type(propensities) == float:
            for sym in self.species:
                propensities *= sp.Pow(sym, reactants.coeff(sym))

        # If propensites in non factorised form, factorise it and generate a dictionary
        if(isinstance(propensities, sp.Expr)):
            after_factor = sp.factor_list(propensities)     # !!! Works only for polynomials right now !!!

            num_factors = len(after_factor[1])
            coefficient = after_factor[0]**(1.0/num_factors)
            
            propensities = {}

            for i in range(num_factors):
                factor = sp.Pow(after_factor[1][i][0],after_factor[1][i][1])
                elements = list(factor.atoms(sp.Symbol))

                if(len(elements) != 1):
                    print("ERROR: Propensity non factorizable")
                    sys.exit()

                propensities[elements[0]] = factor * coefficient

        # Using the dictionary, generate the lambda functions to append the reactions
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



class Partitioning:

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

    def set_initial_condition(self, polynomials_dict):
        
        polynomials = []
        for sym in self.model.species:
            for key, value in list(polynomials_dict.items()):
                if(key == sym):
                    polynomials.append(sp.lambdify(sym,value))

        for Q in self.initial_conditions.Q:
            Q[0, 0, 0] = 1.0

        species_idx = 0
        for node in range(self.tree.n_external_nodes):
            vec_index = np.zeros(self.initial_conditions.external_nodes[node].grid.d())
            for i in range(self.initial_conditions.external_nodes[node].grid.dx()):
                self.initial_conditions.X[node][i, :] = 1
                for j in range(self.initial_conditions.external_nodes[node].grid.d()):
                    self.initial_conditions.X[node][i, :] *= polynomials[species_idx + j](vec_index[j])
                incrVecIndex(vec_index, self.initial_conditions.external_nodes[node].grid.n, self.initial_conditions.external_nodes[node].grid.d())
            species_idx += len(vec_index)
    

def run(partitioning, output, snapshot, tau, tfinal, substeps, method = "implicit_Euler"):
    partitioning.tree.write()
    snap = int(np.floor((tfinal/tau)/snapshot))
    if method == "implicit_Euler":
        m = 'i'
    elif method == "explicit_Euler":
        m = 'e'
    elif method == "Crank-Nicolson":
        m = 'c'
    elif method == "RK4":
        m = 'r'
    else:
        print("Possible inputs for method: implicit_Euler, explicit_Euler, Crank-Nicolson, RK4")
    cmd = f'bin/hierarchical-cme -i input/input.nc -o {output} -s {snap} -t {tau} -f {tfinal} -n {substeps} -m {m}'
    os.system(cmd)