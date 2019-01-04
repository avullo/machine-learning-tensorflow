"""
Define configuration data, utilities and classes to support
the generation of a dataset and the learning problem. 
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

PDB_DIR = '' # where to keep downloaded PDB files
DSSP_DIR = '' # where to keep generated DSSP outputs

DSSP = 'mkdssp' # dssp executable

class SSConfig(object):
    """
    Configs for the Secondary Structure prediction problem.
    """

    pass

class SAConfig(object):
    """
    Configs for the Solvent Accessibility prediction problem.
    """

    pass
