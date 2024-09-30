import unittest

import numpy as np
import sympy as sp

from src.generator import Model, Partitioning, species


class ModelTestCase(unittest.TestCase):
    def setUp(self):
        self.A, self.B = sp.symbols("A,B")
        self.model = Model((self.A, self.B))

    def test_species(self):
        res = species("A")

        self.assertEqual(self.A, res)

    def test_add_reaction(self):
        res_prop = 6

        reactants = 2 * self.A
        products = 3 * self.B
        propensities = 2

        self.model.add_reaction(reactants, products, propensities)
        self.model.generate_reaction_system()

        self.assertEqual(
            self.model.reaction_system.reactions[0].propensity[0](3), res_prop
        )


class PartitioningTestCase(unittest.TestCase):
    def setUp(self):
        self.A, self.B = sp.symbols("A,B")
        self.model = Model((self.A, self.B))

        self.model.add_reaction(2 * self.A, 3 * self.B, 2)

        self.model.generate_reaction_system()

        r = np.array([1])
        p0 = "(A)(B)"
        self.partitioning = Partitioning(p0, r, self.model)

        n = np.array([2, 2])
        d = n.size
        binsize = np.ones(d, dtype=int)
        liml = np.zeros(d)
        self.partitioning.add_grid_params(n, binsize, liml)

        self.partitioning.generate_tree()

        n_basisfunctions = np.ones(r.size, dtype="int")
        self.partitioning.generate_initial_condition(n_basisfunctions)

    def test_set_initial_condition(self):
        res = [[[1]]]

        polynomials_dict = {self.A: self.A * (self.A - 1), self.B: self.B * 2}
        self.partitioning.set_initial_condition(polynomials_dict)

        self.assertEqual(self.partitioning.initial_conditions.Q[0], res)


if __name__ == "__main__":
    unittest.main()
