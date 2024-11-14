import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import scripts.boolean_helper
from scripts.grid_class import GridParms
from scripts.tree_class import Tree, plotReactionGraph

reaction_system = scripts.boolean_helper.convertRulesToReactions("scripts/models/boolean_rulefiles/mTor.hpp")

plt.style.use("./scripts/notebooks/custom_style_boolean.mplstyle")

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

d = 22
n = 2 * np.ones(d, dtype=int)
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

tree = Tree("(0 1 2 3 4 5 6 7 8 9 10)(11 12 13 14 15 16 17 18 19 20 21)", grid)
r_out = np.ones(tree.n_internal_nodes, dtype="int") * 5
tree.initialize(reaction_system, r_out)

plotReactionGraph(tree.G, "plots/mTor_graph.pdf", mode="uniform", color_hex=colors[0], fontcolor_hex="#FFFFFF", fontsize=9, prog="neato", seed=7)