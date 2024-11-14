from distutils.core import setup

setup(name='hierarchical-cme-scripts',
      version='1.0',
      description='Scripts for generating input files and plotting output files',
      packages = ['scripts.models'],
      py_modules = ['scripts.boolean_helper',
                    'scripts.id_class',
                    'scripts.initial_condition_class',
                    'scripts.index_functions',
                    'scripts.grid_class',
                    'scripts.notebooks.output_helper',
                    'scripts.reaction_class',
                    'scripts.tree_class']
      )