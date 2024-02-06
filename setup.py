from distutils.core import setup

setup(name='kinetic-cme-scripts',
      version='1.0',
      description='Scripts for generating input files and plotting output files',
      py_modules = ['scripts.hierarchical.models.bax',
                    'scripts.hierarchical.models.lambda_phage',
                    'scripts.hierarchical.models.toggle_switch',
                    'scripts.hierarchical.id_class',
                    'scripts.hierarchical.initial_condition_class',
                    'scripts.hierarchical.grid_class',
                    'scripts.hierarchical.reaction_class',
                    'scripts.hierarchical.tree_class']
      )
