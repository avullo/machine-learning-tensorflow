"""
Define configuration data, utilities and classes to support the 
generation of a dataset and the definition of the learning problem. 
"""

from __future__ import absolute_import

import os

# package root dir
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# PDB
PDB_URL = 'https://files.rcsb.org/download/'
PDB_DIR = os.path.join( ROOT_DIR, 'data', 'pdb' ) # downloaded PDB files location

# DSSP
DSSP = 'mkdssp' # dssp executable
DSSP_DIR = os.path.join( ROOT_DIR, 'data', 'dssp' ) # location of generated DSSP outputs

class SSConfig(object):
    """
    Configs for the Secondary Structure data generation and prediction problem.
    """

    pass

class SAConfig(object):
    """
    Configs for the Solvent Accessibility data generation and prediction problem.
    """

    pass
