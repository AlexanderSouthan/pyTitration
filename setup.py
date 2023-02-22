from setuptools import setup, find_packages

setup(
    name='pyTitration',
    version='0.0.1',
    packages=find_packages(where='src'),
    install_requires=['numpy', 'pandas', 'matplotlib', 'sympy', 'pip', 'little_helpers']
)
