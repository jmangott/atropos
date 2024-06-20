import numpy as np
import unittest

from scripts.models.bax import reaction_system as bax_model
from scripts.models.lambda_phage import reaction_system as lp_model
from scripts.tree_class import Tree
from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition

class BaxTestCase(unittest.TestCase):
    def setUp(self):
        d = 11
        n = np.arange(1, d+1, dtype=int)
        binsize = np.ones(d)
        liml = np.zeros(d)

        self.partition_str = "(0 1 2)((((3 6)(4 7))(5 8))(9 10))"
        self.r_out = np.array([4, 5, 6, 7])
        self.grid = GridParms(n, binsize, liml)
        self.n_reactions = bax_model.size()

    def test_r_out1(self):
        self.r_out = np.array([4, 5, 6])
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)
            bax_tree.initialize(bax_model, self.r_out)

    def test_r_out2(self):
        self.r_out = np.array([4, 5, 6, 7, 8])
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)
            bax_tree.initialize(bax_model, self.r_out)

    def test_partition_str1(self):
        self.partition_str = "(1 11 2)((((3 6)(4 7))(5 8))(9 10))"
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)

    def test_partition_str2(self):
        self.partition_str = "(1 11 2)(((3 6)(4 7))(5 8))(9 10))"
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)

    def test_partition_str3(self):
        self.partition_str = "(0 1 2)((((3 6)(4 7))(3 8))(9 10))"
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)

    def test_partition_str4(self):
        self.partition_str = "(0 1)((((2 6)(5 7))(3 8))(9 10))"
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)

    def test_reaction_model(self):
        with self.assertRaises(Exception):
            bax_tree = Tree(self.partition_str, self.grid)
            bax_tree.initialize(lp_model, self.r_out)

    def test_bax_tree_partition(self):
        bax_tree = Tree(self.partition_str, self.grid)
        bax_tree.initialize(bax_model, self.r_out)
        
        self.assertEqual(bax_tree.root.grid.dx(), np.prod(self.grid.n))

        self.assertEqual(bax_tree.root.child[0].grid.dx(), np.prod(self.grid.n[:3]))

        self.assertEqual(bax_tree.root.child[1].grid.dx(), np.prod(self.grid.n[3:]))

        self.assertEqual(bax_tree.root.child[1].child[0].grid.dx(), np.prod(self.grid.n[3:9]))

        self.assertEqual(bax_tree.root.child[1].child[1].grid.dx(), np.prod(self.grid.n[9:]))

        self.assertEqual(bax_tree.root.child[1].child[0].child[0].grid.dx(), np.prod(self.grid.n[[3, 6, 4, 7]]))
        self.assertEqual(bax_tree.root.child[1].child[0].child[1].grid.dx(), np.prod(self.grid.n[[5, 8]]))

        self.assertEqual(bax_tree.root.child[1].child[0].child[0].child[0].grid.dx(), np.prod(self.grid.n[[3, 6]]))

        self.assertEqual(bax_tree.root.child[1].child[0].child[0].child[1].grid.dx(), np.prod(self.grid.n[[4, 7]]))
        
        self.assertEqual(bax_tree.root.child[0].child[0], None)

        propensity = [np.array([1.])] * self.n_reactions

        propensity[5] = bax_model.reactions[5].propensity[3](np.arange(self.grid.n[3]))

        propensity[6] = bax_model.reactions[6].propensity[3](np.arange(self.grid.n[3]))

        propensity[10] = bax_model.reactions[10].propensity[3](np.arange(self.grid.n[3]))

        propensity[11] = bax_model.reactions[11].propensity[6](np.arange(self.grid.n[6]))

        propensity[12] = bax_model.reactions[12].propensity[6](np.arange(self.grid.n[6]))

        for i, prop in enumerate(bax_tree.root.child[1].child[0].child[0].child[0].propensity):
            self.assertTrue(np.all(prop == propensity[i]))

        self.dep = np.zeros((4, self.n_reactions), dtype="bool")

        self.dep[0, 5] = True
        self.dep[0, 6] = True
        self.dep[2, 7] = True
        self.dep[2, 8] = True
        self.dep[0, 10] = True
        self.dep[1, 11] = True
        self.dep[1, 12] = True
        self.dep[2, 13] = True
        self.dep[3, 14] = True
        self.dep[3, 15] = True

        self.assertTrue(np.all(bax_tree.root.child[1].child[0].child[0].grid.dep == self.dep))

class LambdaPhageTestCase(unittest.TestCase):
    def setUp(self):
        d = 5
        n = np.array([2, 1, 3, 4, 5])
        binsize = np.ones(d)
        liml = np.zeros(d)
        self.grid = GridParms(n, binsize, liml)
        self.n_reactions = bax_model.size()
        self.partition_str = "((4)(0 1))(2 3)"
        self.r_out = np.array([4, 5])

    def test_lp_partition(self):
        lp_tree = Tree(self.partition_str, self.grid)
        lp_tree.initialize(lp_model, self.r_out)

        self.assertTrue(np.all(lp_tree.root.grid.n == self.grid.n[[4, 0, 1, 2, 3]]))

        self.assertTrue(np.all(lp_tree.root.child[0].grid.n == self.grid.n[[4, 0, 1]]))

        self.assertTrue(np.all(lp_tree.root.child[0].child[0].grid.n == self.grid.n[4]))

        self.assertTrue(np.all(lp_tree.root.child[0].child[1].grid.n == self.grid.n[[0, 1]]))

        self.assertTrue(np.all(lp_tree.root.child[1].grid.n == self.grid.n[2:4]))

        propensity00 = [np.array([1.])] * self.n_reactions

        propensity00[1] = lp_model.reactions[1].propensity[4](np.arange(self.grid.n[4]))

        propensity00[9] = lp_model.reactions[9].propensity[4](np.arange(self.grid.n[4]))

        for i, prop00 in enumerate(lp_tree.root.child[0].child[0].propensity):
            self.assertTrue(np.all(prop00 == propensity00[i]))

        propensity01 = [np.array([1.])] * self.n_reactions

        propensity01[0] = lp_model.reactions[0].propensity[1](np.arange(self.grid.n[1]))

        propensity01[1] = lp_model.reactions[1].propensity[0](np.arange(self.grid.n[0]))

        propensity01[2] = lp_model.reactions[2].propensity[1](np.arange(self.grid.n[1]))

        propensity01[5] = lp_model.reactions[5].propensity[0](np.arange(self.grid.n[0]))

        propensity01[6] = lp_model.reactions[6].propensity[1](np.arange(self.grid.n[1]))

        for i, prop01 in enumerate(lp_tree.root.child[0].child[1].propensity):
            self.assertTrue(np.all(prop01 == propensity01[i]))
        
    def test_slice_vec(self):
        lp_tree = Tree(self.partition_str, self.grid)
        with self.assertRaises(Exception):
            lp_tree.calculateObservables(np.ones(10))

    # TODO: write a better test 
    def test_calculate_observables(self):
        lp_tree = Tree(self.partition_str, self.grid)
        lp_tree.initialize(lp_model, self.r_out)
        initial_condition = InitialCondition(lp_tree, np.ones(self.r_out.size, dtype="int"))
        for node in initial_condition.external_nodes:
            node.X = np.ones((node.grid.dx(), 1))

        for Q in initial_condition.Q:
            Q[0, 0, 0] = 1.0

        slice_vec = np.zeros(lp_tree.grid.d(), dtype="int")
        norm = np.prod(lp_tree.root.grid.n)
        sorted = np.argsort(lp_tree.species)
        n = lp_tree.root.grid.n[sorted]

        sliced_distribution, marginal_distribution = lp_tree.calculateObservables(slice_vec)
        for i, n_el in enumerate(n):
            self.assertTrue(np.all(sliced_distribution[i] == np.ones(n_el)))
            self.assertTrue(np.all(marginal_distribution[i] == np.ones(n_el) * norm / n_el))
        
        with self.assertRaises(Exception):
            _, _ = lp_tree.calculateObservables(np.zeros(lp_tree.grid.d() + 1))

if __name__ == "__main__":
    unittest.main()