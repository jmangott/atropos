import numpy as np
from scripts.reaction_class import Reaction, ReactionSystem

d = 34
reactions = []

# reaction 0

nu1_p = np.zeros(d)
nu1_p[1] = 1
reactions.append(Reaction({1: lambda x: 1-x, 0: lambda x: x}, nu1_p))

nu2_p = np.zeros(d)
nu2_p[2] = 1
reactions.append(Reaction({2: lambda x: 1-x, 0: lambda x: x}, nu2_p))

nu3_p = np.zeros(d)
nu3_p[3] = 1
reactions.append(Reaction({3: lambda x: 1-x, 1: lambda x: x}, nu3_p))

nu4_p = np.zeros(d)
nu4_p[4] = 1
reactions.append(Reaction({4: lambda x: 1-x, 2: lambda x: x}, nu4_p))

nu5_p = np.zeros(d)
nu5_p[5] = 1
reactions.append(Reaction({5: lambda x: 1-x, 3: lambda x: x}, nu5_p))

nu6_p = np.zeros(d)
nu6_p[6] = 1
reactions.append(Reaction({6: lambda x: 1-x, 3: lambda x: x}, nu6_p))

nu7_p = np.zeros(d)
nu7_p[7] = 1
reactions.append(Reaction({7: lambda x: 1-x, 4: lambda x: 1-x, 14: lambda x: x}, nu7_p))
reactions.append(Reaction({7: lambda x: 1-x, 4: lambda x: x, 14: lambda x: 1-x}, nu7_p))
reactions.append(Reaction({7: lambda x: 1-x, 4: lambda x: x, 14: lambda x: x}, nu7_p))

nu8_p = np.zeros(d)
nu8_p[8] = 1
reactions.append(Reaction({8: lambda x: 1-x, 7: lambda x: x}, nu8_p))

nu9_p = np.zeros(d)
nu9_p[9] = 1
reactions.append(Reaction({9: lambda x: 1-x, 5: lambda x: 1-x, 4: lambda x: x}, nu9_p))
reactions.append(Reaction({9: lambda x: 1-x, 5: lambda x: x, 4: lambda x: 1-x}, nu9_p))
reactions.append(Reaction({9: lambda x: 1-x, 5: lambda x: x, 4: lambda x: x}, nu9_p))

nu10_p = np.zeros(d)
nu10_p[10] = 1
reactions.append(Reaction({10: lambda x: 1-x, 6: lambda x: 1-x, 8: lambda x: x}, nu10_p))
reactions.append(Reaction({10: lambda x: 1-x, 6: lambda x: x, 8: lambda x: 1-x}, nu10_p))
reactions.append(Reaction({10: lambda x: 1-x, 6: lambda x: x, 8: lambda x: x}, nu10_p))

nu11_p = np.zeros(d)
nu11_p[11] = 1
reactions.append(Reaction({11: lambda x: 1-x, 10: lambda x: x}, nu11_p))

nu12_p = np.zeros(d)
nu12_p[12] = 1
reactions.append(Reaction({12: lambda x: 1-x, 6: lambda x: x}, nu12_p))

nu13_p = np.zeros(d)
nu13_p[13] = 1
nu13_m = np.zeros(d)
nu13_m[13] = -1
reactions.append(Reaction({13: lambda x: 1-x, 19: lambda x: 1-x, 9: lambda x: x}, nu13_p))
reactions.append(Reaction({13: lambda x: x, 19: lambda x: x, 9: lambda x: 1-x}, nu13_m))
reactions.append(Reaction({13: lambda x: x, 19: lambda x: x, 9: lambda x: x}, nu13_m))

nu14_p = np.zeros(d)
nu14_p[14] = 1
reactions.append(Reaction({14: lambda x: 1-x, 13: lambda x: x}, nu14_p))

nu15_p = np.zeros(d)
nu15_p[15] = 1
reactions.append(Reaction({15: lambda x: 1-x, 10: lambda x: 1-x, 26: lambda x: x}, nu15_p))
reactions.append(Reaction({15: lambda x: 1-x, 10: lambda x: x, 26: lambda x: 1-x}, nu15_p))
reactions.append(Reaction({15: lambda x: 1-x, 10: lambda x: x, 26: lambda x: x}, nu15_p))

# reaction 16

nu17_p = np.zeros(d)
nu17_p[17] = 1
nu17_m = np.zeros(d)
nu17_m[17] = -1

reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: 1-x, 10: lambda x: 1-x, 12: lambda x: x}, nu17_p))
reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: 1-x, 10: lambda x: x, 12: lambda x: 1-x}, nu17_p))
reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: 1-x, 10: lambda x: x, 12: lambda x: x}, nu17_p))

reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: x, 10: lambda x: 1-x, 12: lambda x: 1-x}, nu17_p))
reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: x, 10: lambda x: 1-x, 12: lambda x: x}, nu17_p))
reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: x, 10: lambda x: x, 12: lambda x: 1-x}, nu17_p))
reactions.append(Reaction({17: lambda x: 1-x, 21: lambda x: 1-x, 14: lambda x: x, 10: lambda x: x, 12: lambda x: x}, nu17_p))

reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: 1-x, 10: lambda x: 1-x, 12: lambda x: 1-x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: 1-x, 10: lambda x: 1-x, 12: lambda x: x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: 1-x, 10: lambda x: x, 12: lambda x: 1-x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: 1-x, 10: lambda x: x, 12: lambda x: x}, nu17_m))

reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: x, 10: lambda x: 1-x, 12: lambda x: 1-x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: x, 10: lambda x: 1-x, 12: lambda x: x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: x, 10: lambda x: x, 12: lambda x: 1-x}, nu17_m))
reactions.append(Reaction({17: lambda x: x, 21: lambda x: x, 14: lambda x: x, 10: lambda x: x, 12: lambda x: x}, nu17_m))

nu18_p = np.zeros(d)
nu18_p[18] = 1
nu18_m = np.zeros(d)
nu18_m[18] = -1

reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_p))
reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_p))
reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: x}, nu18_p))

reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_p))
reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_p))
reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_p))
reactions.append(Reaction({18: lambda x: 1-x, 16: lambda x: 1-x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: x}, nu18_p))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: 1-x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: 1-x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: 1-x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: 1-x, 26: lambda x: x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: 1-x}, nu18_m))
reactions.append(Reaction({18: lambda x: x, 16: lambda x: x, 28: lambda x: x, 11: lambda x: x, 15: lambda x: x, 26: lambda x: x}, nu18_m))

nu19_p = np.zeros(d)
nu19_p[19] = 1
reactions.append(Reaction({19: lambda x: 1-x, 25: lambda x: x}, nu19_p))

nu20_p = np.zeros(d)
nu20_p[20] = 1
nu20_m = np.zeros(d)
nu20_m[20] = -1

reactions.append(Reaction({20: lambda x: 1-x, 27: lambda x: 1-x, 14: lambda x: 1-x, 25: lambda x: x}, nu20_p))
reactions.append(Reaction({20: lambda x: 1-x, 27: lambda x: 1-x, 14: lambda x: x, 25: lambda x: 1-x}, nu20_p))
reactions.append(Reaction({20: lambda x: 1-x, 27: lambda x: 1-x, 14: lambda x: x, 25: lambda x: x}, nu20_p))

reactions.append(Reaction({20: lambda x: x, 27: lambda x: x, 14: lambda x: 1-x, 25: lambda x: 1-x}, nu20_m))
reactions.append(Reaction({20: lambda x: x, 27: lambda x: x, 14: lambda x: 1-x, 25: lambda x: x}, nu20_m))
reactions.append(Reaction({20: lambda x: x, 27: lambda x: x, 14: lambda x: x, 25: lambda x: 1-x}, nu20_m))
reactions.append(Reaction({20: lambda x: x, 27: lambda x: x, 14: lambda x: x, 25: lambda x: x}, nu20_m))

nu21_p = np.zeros(d)
nu21_p[21] = 1
reactions.append(Reaction({21: lambda x: 1-x, 26: lambda x: x}, nu21_p))

nu22_p = np.zeros(d)
nu22_p[22] = 1
nu22_m = np.zeros(d)
nu22_m[22] = -1
reactions.append(Reaction({22: lambda x: 1-x, 24: lambda x: 1-x, 15: lambda x: x}, nu22_p))
reactions.append(Reaction({22: lambda x: x, 24: lambda x: x, 15: lambda x: 1-x}, nu22_m))
reactions.append(Reaction({22: lambda x: x, 24: lambda x: x, 15: lambda x: x}, nu22_m))

nu23_p = np.zeros(d)
nu23_p[23] = 1
nu23_m = np.zeros(d)
nu23_m[23] = -1
reactions.append(Reaction({23: lambda x: 1-x, 17: lambda x: 1-x, 26: lambda x: x}, nu23_p))
reactions.append(Reaction({23: lambda x: x, 17: lambda x: x, 26: lambda x: 1-x}, nu23_m))
reactions.append(Reaction({23: lambda x: x, 17: lambda x: x, 26: lambda x: x}, nu23_m))

nu24_m = np.zeros(d)
nu24_m[24] = -1
reactions.append(Reaction({24: lambda x: x, 18: lambda x: 1-x, 31: lambda x: x}, nu24_m))
reactions.append(Reaction({24: lambda x: x, 18: lambda x: x, 31: lambda x: 1-x}, nu24_m))
reactions.append(Reaction({24: lambda x: x, 18: lambda x: x, 31: lambda x: x}, nu24_m))

nu25_m = np.zeros(d)
nu25_m[25] = -1
reactions.append(Reaction({25: lambda x: x, 20: lambda x: x}, nu25_m))

nu26_m = np.zeros(d)
nu26_m[26] = -1
reactions.append(Reaction({26: lambda x: x, 23: lambda x: x}, nu26_m))

nu27_p = np.zeros(d)
nu27_p[27] = 1
reactions.append(Reaction({27: lambda x: 1-x, 22: lambda x: x}, nu27_p))

nu28_p = np.zeros(d)
nu28_p[28] = 1
reactions.append(Reaction({28: lambda x: 1-x, 25: lambda x: x}, nu28_p))

nu29_p = np.zeros(d)
nu29_p[29] = 1
reactions.append(Reaction({29: lambda x: 1-x, 25: lambda x: x}, nu29_p))

nu30_p = np.zeros(d)
nu30_p[30] = 1
nu30_m = np.zeros(d)
nu30_m[30] = -1
reactions.append(Reaction({30: lambda x: 1-x, 25: lambda x: 1-x, 26: lambda x: x}, nu30_p))
reactions.append(Reaction({30: lambda x: x, 25: lambda x: x, 26: lambda x: 1-x}, nu30_m))
reactions.append(Reaction({30: lambda x: x, 25: lambda x: x, 26: lambda x: x}, nu30_m))

nu31_p = np.zeros(d)
nu31_p[31] = 1
nu31_m = np.zeros(d)
nu31_m[31] = -1
reactions.append(Reaction({31: lambda x: 1-x, 28: lambda x: 1-x, 22: lambda x: x}, nu31_p))
reactions.append(Reaction({31: lambda x: x, 28: lambda x: x, 22: lambda x: 1-x}, nu31_m))
reactions.append(Reaction({31: lambda x: x, 28: lambda x: x, 22: lambda x: x}, nu31_m))

nu32_p = np.zeros(d)
nu32_p[32] = 1
nu32_m = np.zeros(d)
nu32_m[32] = -1
reactions.append(Reaction({32: lambda x: 1-x, 30: lambda x: 1-x, 29: lambda x: x}, nu32_p))
reactions.append(Reaction({32: lambda x: x, 30: lambda x: x, 29: lambda x: 1-x}, nu32_m))
reactions.append(Reaction({32: lambda x: x, 30: lambda x: x, 29: lambda x: x}, nu32_m))

nu33_m = np.zeros(d)
nu33_m[33] = -1
reactions.append(Reaction({33: lambda x: x, 31: lambda x: 1-x}, nu33_m))

species_names = ["HMGB1","TLR24","RAGE","MYD88","RAS","RAC1","IRAKs","RAF","MEK","PI3K","ERK","AP1","TAB1","PIP3","AKT","Myc","INK4a",
                 "IKK","CyclinD","PTEN","MDM2","A20","E2F","IkB","RB","P53","NFkB","ARF","P21","BAX","BclXL","CyclinE","Apoptosis","Proliferate"]

reaction_system = ReactionSystem(reactions, species_names)

if __name__ == "__main__":
    import collections
    import itertools
    import networkx as nx
    import matplotlib.pyplot as plt
    reactants = [reaction.propensity.keys() for reaction in reaction_system.reactions]
    combinations = [comb for reactant in reactants for comb in itertools.combinations(reactant, 2)]
    counter = collections.Counter(combinations)
    edges = counter.keys()
    weights = counter.values()
    edges_weights = [(e[0], e[1], {"weight": w}) for e, w in zip(edges, weights)]
    G = nx.Graph(edges_weights)
    nx.set_node_attributes(G, {i: name for i, name in enumerate(reaction_system.species_names)}, "labels")

    # Plot graph
    widths = np.fromiter(nx.get_edge_attributes(G, 'weight').values(), dtype=float)
    pos = nx.spring_layout(G)

    plt.figure(figsize=(12,8))
    nx.draw_networkx_nodes(G, pos, node_size=750)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=np.log10(widths*10), alpha=0.6)
    nx.draw_networkx_labels(G, pos, nx.get_node_attributes(G, "labels"), font_size=8)
    plt.show()