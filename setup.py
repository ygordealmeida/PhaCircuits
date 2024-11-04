from setuptools import setup, find_packages

setup(
    name="PhaCircuits",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.24.1",
        "schemdraw>=0.19",
    ],
    author="@ygordealmeida",
    description="Circuits simulation in phasorial domain",
)
