from distutils.core import setup

setup(
    name="hierarchical-cme-scripts",
    version="1.0",
    description="Scripts for generating input files and plotting output files",
    packages=[
        "examples",
        "examples.boolean",
        "examples.kinetic",
        "examples.models.kinetic",
        "src",
    ],
)
