from setuptools import setup, find_packages

setup(
    name="dedaLES",
    packages=find_packages(),
    use_scm_version=True,
    setup_requires=['setuptools_scm'])

