from unittest import TestCase

from protstruct.dataset import get_pdb_file
from protstruct.config import PDB_DIR

import os
import tensorflow as tf

class TestGetPDBFile(TestCase):
    """
    Test routine for downloading PDB files
    """

    def test_uncompressed_file(self):
        filename = get_pdb_file( '12as' )

        self.assertEqual( filename, os.path.join( PDB_DIR, '12as.pdb' ) )
        self.assertTrue( tf.gfile.Exists( filename ) )

        tf.gfile.Remove( filename )

    def test_compressed_file(self):
        filename = get_pdb_file( '16vp', compressed=True )

        self.assertEqual( filename, os.path.join( PDB_DIR, '16vp.pdb.gz' ) )
        self.assertTrue( tf.gfile.Exists( filename ) )

        tf.gfile.Remove( filename )
