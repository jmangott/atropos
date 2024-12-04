import numpy as np
import unittest

import scripts.boolean_helper
from scripts.tree_class import Tree
from scripts.grid_class import GridParms

class PancreaticCancerTestCase(unittest.TestCase):
    def setUp(self):
        d = 34
        n = 2 * np.ones(d, dtype=int)
        binsize = np.ones(d, dtype=int)
        liml = np.zeros(d)

        self.reaction_system = scripts.boolean_helper.convertRulesToReactions("scripts/models/boolean_rulefiles/pancreatic_cancer.hpp")
        self.grid = GridParms(n, binsize, liml)

    def test_one_level(self):
        p_best = "(0 1 2 3 4 5 6 7 8 9 10 11 12 17 21 23 26)(13 14 15 16 18 19 20 22 24 25 27 28 29 30 31 32 33)"
        p_worst = "(0 1 2 11 15 16 18 19 20 21 23 25 26 28 29 30 32)(3 4 5 6 7 8 9 10 12 13 14 17 22 24 27 31 33)"
        p_reasonable = "(0 1 2 3 4 5 7 9 13 14 19 20 25 27 29 30 32)(6 8 10 11 12 15 16 17 18 21 22 23 24 26 28 31 33)"
        p_literature = "(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)(17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33)"

        partition_strings = [p_best, p_worst, p_reasonable, p_literature]

        reference_entropies = [1.26884395, 6.39031953, 2.3125, 2.00519536]

        for reference_entropy, partition_str in zip(reference_entropies, partition_strings):
            tree = Tree(partition_str, self.grid)
            r_out = np.ones(tree.n_internal_nodes, dtype="int") * 1
            tree.initialize(self.reaction_system, r_out)
            entropy = tree.calculateEntropy(tree.root)

            self.assertAlmostEqual(entropy, reference_entropy)

    def test_hierarchical(self):
        p_best = "((0 1 2 3 4 5 7 8 9)(6 10 11 12 17 21 23 26))((13 14 19 20 25 27 29 30 32)(15 16 18 22 24 28 31 33))"
        p_worst = "((0 1 19 20 21 23 25 26 30)(2 11 15 16 18 28 29 32))((3 4 5 6 8 9 10 13)(7 12 14 17 22 24 27 31 33))"
        p_reasonable = "((0 1 2 3 4 5 7 9)(13 14 19 20 25 27 29 30 32))((6 8 10 12 17 21 23 26)(11 15 16 18 22 24 28 31 33))"

        partition_strings = [p_best, p_worst, p_reasonable]

        reference_entropies = [3.01884395, 11.04656953, 3.74938906]

        for reference_entropy, partition_str in zip(reference_entropies, partition_strings):
            tree = Tree(partition_str, self.grid)
            r_out = np.ones(tree.n_internal_nodes, dtype="int") * 1
            tree.initialize(self.reaction_system, r_out)
            entropy = (
                tree.calculateEntropy(tree.root)
                + tree.calculateEntropy(tree.root.child[0])
                + tree.calculateEntropy(tree.root.child[1])
            )

            self.assertAlmostEqual(entropy, reference_entropy)

class ApoptosisTestCase(unittest.TestCase):
    def setUp(self):
        d = 41
        n = 2 * np.ones(d, dtype=int)
        binsize = np.ones(d, dtype=int)
        liml = np.zeros(d)

        self.reaction_system = scripts.boolean_helper.convertRulesToReactions(
            "scripts/models/boolean_rulefiles/apoptosis.hpp"
        )
        self.grid = GridParms(n, binsize, liml)

    def test_hierarchical(self):
        p_best = "((0 2 3 4 5 6 7 8 9 12 20)(21 25 26 27 28 29 30 31 39 40))((1 10 11 14 15 22 23 32 33 38)(13 16 17 18 19 24 34 35 36 37))"
        p_worst = "((0 2 4 12 13 19 24 32 35 37)(3 5 7 8 11 22 23 33 34 38))((1 6 9 10 14 20 21 25 26 36)(15 16 17 18 27 28 29 30 31 39 40))"
        p_reasonable = "((0 2 3 4 5 6 7 8 12 20)(1 9 10 11 14 15 16 17 22 23 37))((13 18 19 24 32 33 34 35 36 38)(21 25 26 27 28 29 30 31 39 40))"

        partition_strings = [p_best, p_worst, p_reasonable]

        reference_entropies = [5.97396601, 18.99050742, 7.24877812]

        for reference_entropy, partition_str in zip(
            reference_entropies, partition_strings
        ):
            tree = Tree(partition_str, self.grid)
            r_out = np.ones(tree.n_internal_nodes, dtype="int") * 1
            tree.initialize(self.reaction_system, r_out)
            entropy = (
                tree.calculateEntropy(tree.root)
                + tree.calculateEntropy(tree.root.child[0])
                + tree.calculateEntropy(tree.root.child[1])
            )

            self.assertAlmostEqual(entropy, reference_entropy)

if __name__ == "__main__":
    unittest.main()