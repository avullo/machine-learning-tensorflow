"""
Utilities to download, preprocess datasets.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from utils import download_file

import tensorflow as tf
import config

def get_pdb_file( pdb_file, compressed = True ):
    pdb_url, pdb_dir = config.PDB_URL, config.PDB_DIR

    if not tf.gfile.Exists( pdb_dir ):
        tf.gfile.MakeDirs( pdb_dir )

    pdb_file = os.path.join( pdb_url, pdb_file )
    if compressed:
        pdb_file += '.gz'

    # TODO: test creation of remote path and if method works
    download_file( pdb_url + pdb_file, pdb_dir )

