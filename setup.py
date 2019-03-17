import setuptools
from setuptools import setup

# Grab long_description from README
with open('README.md') as f:
    long_description = f.read()

setup(
    name = "dedaLES",
    author = "Keaton Burns, Gregory L. Wagner",
    author_email = "keaton.burns@gmail.com, wagner.greg@gmail.com",
    description = "A package for performing Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES) of fluid flows using the Dedalus framework.",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/glwagner/dedaLES",
    classifiers = ["Programming Language :: Python :: 3"],
    packages = setuptools.find_packages(),
    use_scm_version = True,
    setup_requires = ['setuptools_scm'])
