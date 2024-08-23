"""
Contains the `Tree` class, which stores the low-rank approximation of an initial probability distribution as a binary tree according to a prescribed partition.
"""
import collections
import copy
from datatree import DataTree
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import os
import regex
from scipy.special import xlogy
import xarray as xr

from scripts.grid_class import GridParms
from scripts.id_class import Id
from scripts.index_functions import incrVecIndex, vecIndexToCombIndex, vecIndexToState, tensorUnfold
from scripts.reaction_class import ReactionSystem

class Node:
    def __init__(self, _id: Id, _grid: GridParms):
        self.child = [None] * 2
        self.id = _id
        self.grid = _grid

class InternalNode(Node):
    def __init__(self, _parent: 'InternalNode', _id: Id, _grid: GridParms):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.Q = np.zeros((0, 0, 0))

    def rankIn(self):
        return self.Q.shape[-1]

    def rankOut(self):
        return self.Q.shape[0]

class ExternalNode(Node):
    def __init__(self, _parent: InternalNode, _id: Id, _grid: GridParms):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.X = np.zeros((0, 0))
        self.propensity = []

    def rankIn(self):
        return self.X.shape[-1]


"""
The `tree` class stores the initial condition for the hierarchical DLR approximation of the chemical master equation. Note that the input for the grid parameters (`_grid`) must follow the same ordering convention for the species as the reaction system (`_reaction_system`), when initializing the tree via the `initialize` method, but the partition string for dividing the reaction networks allows also permutation.
"""
class Tree:
    # TODO: Move the `_grid` dependency to `initialize` method
    def __init__(self,
                 _partition_str: str,
                 _grid: GridParms
                ):

        # Syntactic check (test whether `_partition_str` is a valid input string)
        if not regex.fullmatch(r"\((?:\d+(?:\s\d)*|(?R))+\)\((?:\d+(?:\s\d)*|(?R))+\)", _partition_str):
            raise Exception("Not a valid `_partition_str`")

        p = self.__removeBrackets(_partition_str)
        self.species = np.copy(p)  # Create a deep copy of `p` before sorting it
        p.sort()
        p_diff = np.diff(p)

        # Semantic checks
        if p[0] != 0:
            raise Exception("The smallest value of `_partition_str` has to be 0")
        if np.any(p_diff != 1):
            raise Exception("Not all species covered by `_partition_str`")

        # Test whether the dimension of `_partition_str` and `_grid.d` match
        if p.size != _grid.d():
            raise Exception(
                "Dimensions of `_partition_str` and `_grid.d` do not match")

        self.n_internal_nodes = 1  # 1 is for the root node
        self.n_external_nodes = 0

        self.internal_nodes = {}
        self.external_nodes = {}

        self.reaction_system = None
        self.G = None
        self.species_names = None
        self.partition_str = _partition_str
        self.grid = copy.deepcopy(_grid)
        root_id = Id("")
        self.root = InternalNode(None, root_id, self.grid.permute(self.species))
        self.internal_nodes[root_id] = self.root
        self.__build(self.root, self.partition_str)

    @staticmethod
    def __parsingHelper(input_str):
        if (input_str == "("):
            sigma = 1
        elif (input_str == ")"):
            sigma = -1
        else:
            sigma = 0
        return sigma

    @staticmethod
    def __removeBrackets(input_str):
        n = input_str.replace("(", " ").replace(")", " ").split()
        n = np.array([int(ele) for ele in n])
        return n

    def __build(self, node: InternalNode, partition_str):
        sigma = 0
        i = 0
        for i, ele in enumerate(partition_str):
            sigma += self.__parsingHelper(ele)
            if sigma == 0:
                break

        partition_str0 = partition_str[1:i]
        partition_str1 = partition_str[i+2:-1]
        p0_size = self.__removeBrackets(partition_str0).size

        grid0 = GridParms(
            node.grid.n[:p0_size], node.grid.binsize[:p0_size], node.grid.liml[:p0_size], node.grid.species[:p0_size])
        grid1 = GridParms(
            node.grid.n[p0_size:], node.grid.binsize[p0_size:], node.grid.liml[p0_size:], node.grid.species[p0_size:])
                
        if (partition_str0[0] == "("):
            new_id = node.id + 0
            node.child[0] = InternalNode(node, new_id, grid0)
            self.__build(node.child[0], partition_str0)
            self.n_internal_nodes += 1
            self.internal_nodes[new_id] = node.child[0]
        else:
            new_id = node.id + 0
            node.child[0] = ExternalNode(node, new_id, grid0)
            self.n_external_nodes += 1
            self.external_nodes[new_id] = node.child[0]

        if (partition_str1[0] == "("):
            new_id = node.id + 1
            node.child[1] = InternalNode(node, new_id, grid1)
            self.__build(node.child[1], partition_str1)
            self.n_internal_nodes += 1
            self.internal_nodes[new_id] = node.child[1]
        else:
            new_id = node.id + 1
            node.child[1] = ExternalNode(node, new_id, grid1)
            self.n_external_nodes += 1
            self.external_nodes[new_id] = node.child[1]

        return
    
    def __initialize(self, node: Node):
        p0_size = node.parent.child[0].grid.d()

        idx = node.parent.child.index(node)
        if idx == 0:
            sl = slice(0,p0_size)
        elif idx == 1:
            sl = slice(p0_size, node.parent.grid.d())

        node.grid = GridParms(
            node.parent.grid.n[sl], node.parent.grid.binsize[sl], node.parent.grid.liml[sl], node.parent.grid.species[sl], node.parent.grid.dep[sl, :], node.parent.grid.nu[sl, :])

        if isinstance(node, InternalNode):
            next_r_out = next(self.__r_out_iter)
            node.Q.resize((next_r_out, next_r_out, node.parent.rankOut()))
            self.__initialize(node.child[0])
            self.__initialize(node.child[1])

        if isinstance(node, ExternalNode):
            node.X.resize((node.grid.dx(), node.parent.rankOut()))
            node.propensity = self.__calculatePropensity(node)

        return

    def initialize(self, reaction_system: ReactionSystem, r_out: npt.NDArray[np.int_]):
        # Test whether the dimension of `_r` is equal to n_internal_nodes
        if r_out.size != self.n_internal_nodes:
            raise Exception(
                "`r_out.size` must be equal to the number of internal nodes")
        
        if self.grid.d() != reaction_system.d():
            raise Exception(
                "`self.grid.d()` must be equal to the number of species in the reaction system")

        self.reaction_system = reaction_system
        self.grid.initialize(reaction_system)
        self.root.grid = self.grid.permute(self.species)
        self.__r_out_iter = iter(r_out)
        next_r_out = next(self.__r_out_iter)
        self.root.Q.resize((next_r_out, next_r_out, 1))
        self.__initialize(self.root.child[0])
        self.__initialize(self.root.child[1])
        self.species_names = self.reaction_system.species_names
        self.G = self.__getReactionGraph()

    def __print(self, node: Node, os: str) -> str:
        os = " ".join([os, str(type(node)), "id:", str(node.id), "n:", str(node.grid), "species:", str(node.grid.species)])
        if isinstance(node, ExternalNode):
            os = " ".join([os, "X.shape:", str(node.X.shape), "\n"])
        elif isinstance(node, InternalNode):
            os = " ".join([os, "Q.shape:", str(node.Q.shape), "\n"])
        if node.child[0]:
            os = self.__print(node.child[0], os)
        if node.child[1]:
            os = self.__print(node.child[1], os)
        return os

    def __str__(self) -> str:
        return self.__print(self.root, "")

    @staticmethod
    def __createDataset(node: Node):
        ds = xr.Dataset(
            {
                "n": (["d"], node.grid.n),
                "binsize": (["d"], node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "species": (["d"], node.grid.species),
                "dep": (["d", "n_reactions"], node.grid.dep),
                "nu": (["d", "n_reactions"], node.grid.nu)
            }
        )
        return ds

    def __write(self, node: Node, parent_dt: DataTree):
        if isinstance(node, ExternalNode):
            ds = self.__createDataset(node)
            ds["X"] = (["n_basisfunctions", "dx"], node.X.T)
            for mu, propensity in enumerate(node.propensity):
                ds["propensity_{}".format(mu)] = (["dx_{}".format(mu)], propensity)
            DataTree(name=str(node.id),
                     parent=parent_dt, data=ds)

        elif isinstance(node, InternalNode):
            ds = self.__createDataset(node)
            ds["Q"] = (["n_basisfunctions", "r_out0", "r_out1"], node.Q.T)
            dt = DataTree(name=str(node.id),
                          parent=parent_dt, data=ds)
            self.__write(node.child[0], dt)
            self.__write(node.child[1], dt)

        return

    def write(self, fname: str = "input/input.nc"):
        if not os.path.exists("input"):
            os.makedirs("input")

        # Undo permutation of grid for root only
        self.grid.permute(self.species)

        ds = self.__createDataset(self.root)
        ds["species_names"] = (["d"], self.species_names)
        ds["Q"] = (["n_basisfunctions", "r_out0", "r_out1"], self.root.Q.T)
        dt = DataTree(name=str(self.root.id), data=ds)
        dt.attrs["partition_str"] = self.partition_str
        self.__write(self.root.child[0], dt)
        self.__write(self.root.child[1], dt)
        dt.to_netcdf(fname, engine="netcdf4")

    def __calculateObservableHelper(self, node: Node, idx_n: int, slice_vec: npt.NDArray[np.int_]):
        if isinstance(node, ExternalNode):
            if idx_n in node.grid.species:
                # Sliced distribution
                partition_idx_n = np.where(idx_n==node.grid.species)[0][0]
                sliced_distribution = np.zeros((node.grid.n[partition_idx_n], node.rankIn()), dtype="float64")
                partition_slice_vec = slice_vec[node.grid.species]
                for i in range(node.grid.n[partition_idx_n]):
                    partition_slice_vec[partition_idx_n] = i
                    sliced_distribution[i, :] = node.X[vecIndexToCombIndex(partition_slice_vec, node.grid.n), :]
                # Marginal distribution
                vec_index = np.zeros(node.grid.d(), dtype=np.int64)
                marginal_distribution = np.zeros((node.grid.n[partition_idx_n], node.rankIn()), dtype="float64")
                for i in range(node.grid.dx()):
                    marginal_distribution[vec_index[partition_idx_n], :] += node.X[i, :]
                    incrVecIndex(vec_index, node.grid.n, node.grid.d())
            else:
                partition_slice_vec = slice_vec[node.grid.species]
                sliced_distribution = node.X[vecIndexToCombIndex(partition_slice_vec, node.grid.n), :]
                marginal_distribution = np.sum(node.X, axis=0)
        elif isinstance(node, InternalNode):
            X0_sliced, X0_marginal = self.__calculateObservableHelper(node.child[0], idx_n, slice_vec)
            X1_sliced, X1_marginal = self.__calculateObservableHelper(node.child[1], idx_n, slice_vec)
            if (X0_sliced.ndim > 1):
                rule = "ij,jkl,k->il"
            elif (X1_sliced.ndim > 1):
                rule = "i,ijk,lj->lk"
            else:
                rule = "i,ijk,j->k"
            sliced_distribution = np.einsum(rule, X0_sliced, node.Q, X1_sliced)
            marginal_distribution = np.einsum(rule, X0_marginal, node.Q, X1_marginal)
        return sliced_distribution, marginal_distribution

    def __calculateObservable(self, idx_n: int, slice_vec: npt.NDArray[np.int_]):
        sliced_distribution, marginal_distribution = self.__calculateObservableHelper(self.root, idx_n, slice_vec)
        return sliced_distribution[:, 0], marginal_distribution[:, 0]

    def calculateObservables(self, slice_vec: npt.NDArray[np.int_]):
        sliced_distributions = {}
        marginal_distributions = {}
        for i in self.species:
            sliced, marginal = self.__calculateObservable(i, slice_vec)
            sliced_distributions[self.species_names[i]] = sliced
            marginal_distributions[self.species_names[i]] = marginal
        return sliced_distributions, marginal_distributions
    
    def __calculateFullDistributionHelper(self, node: Node):
        if isinstance(node, ExternalNode):
            X = node.X
        elif isinstance(node, InternalNode):
            X0 = self.__calculateFullDistributionHelper(node.child[0])
            X1 = self.__calculateFullDistributionHelper(node.child[1])
            X_tensor = np.einsum("lmk,il,jm->ijk", node.Q, X0, X1)
            X = tensorUnfold(X_tensor, 2).T
        return X
    
    def calculateFullDistribution(self):
        """
        This method only works when all species occur in ascending order in the partition string. The full probability distribution is computed, therefore this method should be used only for small system sizes. 
        """
        return self.__calculateFullDistributionHelper(self.root)[:, 0]

    def __getReactionGraph(self):
        combinations = []
        reaction_dependencies = self.__getReactionDependencies()

        for (reactants, products) in reaction_dependencies.keys():
            for reactant in reactants:
                for product in products:
                    if reactant != product:
                        combinations.append((reactant, product))

        # combinations = [comb for reactant in reactants for comb in itertools.combinations(reactant, 2)]
        counter = collections.Counter(combinations)
        edges = counter.keys()
        weights = counter.values()
        edges_weights = [(e[0], e[1], {"weight": w}) for e, w in zip(edges, weights)]

        # total_species = [e.grid.species for e in self.external_nodes]
        external_ids = self.external_nodes.keys()
        attributes = {self.species_names[species]: {"id": id} for id in external_ids for species in self.external_nodes[id].grid.species}

        G = nx.Graph(edges_weights)
        nx.set_node_attributes(G, attributes)
        return G

    def __calculatePropensity(self, node: Node):
        propensity = [None] * self.reaction_system.size()
        for mu, reaction in enumerate(self.reaction_system.reactions):
            n_dep = node.grid.n[node.grid.dep[:, mu]]
            dx_dep = np.prod(n_dep)
            propensity[mu] = np.ones(dx_dep)
            vec_index = np.zeros(n_dep.size)
            reactants = [reactant for reactant in node.grid.species if reactant in reaction.propensity.keys()]
            for i in range(dx_dep):
                for j, reactant in enumerate(reactants):
                    propensity[mu][i] *= reaction.propensity[reactant](vec_index[j])
                incrVecIndex(vec_index, n_dep, n_dep.size)
        return propensity

    def __getReactionDependencies(self):
        """
        This method computes the input/output dependencies of the reaction network as a dictionary.
        """
        reaction_dependencies = collections.defaultdict(list)
        for mu, reaction in enumerate(self.reaction_system.reactions):
            input = tuple([self.species_names[k] for k in reaction.propensity.keys()])
            output = tuple([self.species_names[k] for k in np.nonzero(reaction.nu)[0]])
            reaction_dependencies[(input, output)].append(mu)
        return reaction_dependencies

    # TODO: check if this function also works for general kinetic models
    # TODO: for kinetic models the products may be more than one, 
    # but then the reactions have to be converted to a 'normal' form, where only a single product is present
    def calculateEntropy(self, node: InternalNode):
        reaction_dependencies = self.__getReactionDependencies()
        S = {}

        for key, val in reaction_dependencies.items():
            S[key] = 0

            # take the first reaction (all other reactions have the same `grid.dep` values)
            mu_0 = reaction_dependencies[key][0]

            # mapping species_names <-> species_id for the products
            products = findIndex(self.species_names, key[1])
            
            # TODO: rewrite this function that it is clearly visible that only a single product is allowed
            prod_lies_in_partition_0 = (products[0] in node.child[0].grid.species)
            if prod_lies_in_partition_0:
                grid0 = node.child[0].grid
                grid1 = node.child[1].grid
                propensity0 = self.__calculatePropensity(node.child[0])
                propensity1 = self.__calculatePropensity(node.child[1])
            else: # exchange the partitions
                grid0 = node.child[1].grid
                grid1 = node.child[0].grid
                propensity0 = self.__calculatePropensity(node.child[1])
                propensity1 = self.__calculatePropensity(node.child[0])

            d0 = grid0.d()
            d1 = grid1.d()
            dep0 = grid0.dep[:, mu_0]
            dep1 = grid1.dep[:, mu_0]
            n_dep0 = grid0.n[grid0.dep[:, mu_0]]
            n_dep1 = grid1.n[grid1.dep[:, mu_0]]
            dx_dep0 = np.prod(grid0.n[grid0.dep[:, mu_0]])
            dx_dep1 = np.prod(grid1.n[grid1.dep[:, mu_0]])

            species0 = list(grid0.species)
            species1 = list(grid1.species)

            # obtain species_id for the products in partition 0 and 1
            products0 = findIndex(species0, products)
            products1 = findIndex(species1, products)

            state0 = np.zeros(d0)
            dep_vec_index0 = np.zeros(n_dep0.size, dtype=int)

            for i0 in range(dx_dep0):
                count = collections.defaultdict(int)
                state1 = np.zeros(d1)
                dep_vec_index1 = np.zeros(n_dep1.size, dtype=int)
                for i1 in range(dx_dep1):
                    product_population_number0 = state0[products0]
                    product_population_number1 = state1[products1]
                    weight = 1.0
                    for mu in val: # `val` are all reactions with the same input/output dependencies
                        nu0 = self.reaction_system.reactions[mu].nu[species0]
                        nu1 = self.reaction_system.reactions[mu].nu[species1]
                        propensity = propensity0[mu][i0] * propensity1[mu][i1]
                        if not np.isclose(propensity, 0.0, atol=1e-12):
                            # weight *= propensity # for the kinetic case?
                            product_population_number0 += nu0[products0]
                            product_population_number1 += nu1[products1]

                    count[(tuple(product_population_number0), tuple(product_population_number1))] += weight 
                    incrVecIndex(dep_vec_index1, n_dep1, n_dep1.size)
                    state1[dep1] = dep_vec_index1

                count_values = np.array(list(count.values()))
                sum_count_values = np.sum(count_values)
                if sum_count_values != 0:
                    probabilities = count_values / sum_count_values
                    # NOTE: `xlogy` handles the special case `probability=0`
                    S[key] -= np.sum(xlogy(probabilities, probabilities)) / np.log(2.0)
                incrVecIndex(dep_vec_index0, n_dep0, n_dep0.size)
                state0[dep0] = dep_vec_index0
            S[key] /= dx_dep0

        total_entropy = np.sum(np.array(list(S.values())))

        return total_entropy

def findIndex(array: list, values):
    idx = []
    for v in values:
        try:
            idx.append(array.index(v))
        except ValueError:
            pass
    return idx

def plotReactionGraph(G: nx.Graph, seed=None):
    """
    Helper function for plotting the `nx.Graph` member variable `G` of a `Tree` object.
    """
    widths = np.fromiter(nx.get_edge_attributes(G, 'weight').values(), dtype=float)
    pos = nx.spring_layout(G, seed=seed)
    fig, ax = plt.subplots()

    id = nx.get_node_attributes(G, "id")
    color_id = [int(id[id_key]) for id_key in id]
    nx.draw_networkx(G, pos, with_labels=True, node_size=750, ax=ax, node_color=color_id, cmap="tab20")
    nx.draw_networkx_edges(G, pos, width=np.log10(widths*10), alpha=0.6, ax=ax)
    return fig, ax


if __name__ == "__main__":
    import scripts.models.boolean_pancreatic_cancer as model
    d = 34
    n = 2 * np.ones(d, dtype=int)
    binsize = np.ones(d, dtype=int)
    liml = np.zeros(d)
    grid = GridParms(n, binsize, liml)

    partition_str = '(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)(17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33)'

    tree = Tree(partition_str, grid)
    r_out = np.ones(tree.n_internal_nodes, dtype="int") * 5

    tree.initialize(model.reaction_system, r_out)
    entropy = tree.calculateEntropy(tree.root)
    # for i, s in enumerate(tree.species_names):
    #     print(i, s)
    print(entropy)

    fig, ax = plotReactionGraph(tree.G, seed=42)
    plt.show()