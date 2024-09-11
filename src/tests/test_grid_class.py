import unittest

import numpy as np

from examples.models.kinetic.lambda_phage import reaction_system as lp_model
from src.grid import GridParms

class GridTestCase(unittest.TestCase):
    def setUp(self):
        self.n = np.array([1, 2, 3, 4, 5])
        self.binsize = np.array([1, 1, 1, 1, 1])
        self.liml = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    def test_grid_binsize(self):
        self.binsize = np.array([1, 1, 1, 1])
        with self.assertRaises(Exception):
            grid = GridParms(self.n, self.binsize, self.liml)
 
    def test_grid_n_zero(self):
        self.n = np.array([1, 2, 3, 0, 5])
        with self.assertRaises(Exception):
            grid = GridParms(self.n, self.binsize, self.liml)

    def test_grid_initialize(self):
        grid = GridParms(self.n, self.binsize, self.liml)
        grid.initialize(lp_model)
        
        dep = np.zeros((5, lp_model.size()), dtype="bool")
        dep[1, 0] = True
        dep[0, 1] = True
        dep[4, 1] = True
        dep[1, 2] = True
        dep[2, 3] = True
        dep[2, 4] = True
        dep[0, 5] = True
        dep[1, 6] = True
        dep[2, 7] = True
        dep[3, 8] = True
        dep[4, 9] = True

        nu = np.zeros((5, lp_model.size()))
        for i in range(5):
            nu[i, i] = 1
            nu[i, i+5] =-1

        self.assertTrue(np.all(grid.dep == dep))

        self.assertTrue(np.all(grid.nu == nu))

if __name__ == "__main__":
    unittest.main()