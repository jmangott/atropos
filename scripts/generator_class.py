from scripts.reaction_class import Reaction, ReactionSystem
from scripts.grid_class import GridParms
from scripts.tree_class import Tree
import sympy as sp
import numpy as np

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
                    prop_dict[i] = lambda x: value.subs(sym,x)
        
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

    
        


def run(model, partitioning):
    pass




"""
Define some symbols, to then try out with those symbols
"""
NF, GR, O, H = sp.symbols("NF, GR, O, H")



"""
Try out model
"""
model = Model((NF, GR, O, H))
print(model.species, model.reactions)

model.add_reaction(3*NF + 7*GR, 2*H + O, {NF: NF**2, H: 1/(1 + H**2)})
print(model.reactions)

model.add_reactions([7*H + GR, 2*H],[NF, 3*O + 2*NF],[{H: H}, {H: 3*H}])        #TODO: Some error here, dictionary used in tree generation only works with last species
print(model.reactions)

model.generate_reaction_system()
print(model.reaction_system)



"""
Try out partitioning
"""
r = np.array([5,5])
partitioning = Partitioning('((NF GR)(O))(H)', r, model)
print(partitioning.partition)

n = np.array([16, 11, 11, 11])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
partitioning.add_grid_params(n, binsize, liml)
print(partitioning.grid)

partitioning.generate_tree()
print(partitioning.tree)

