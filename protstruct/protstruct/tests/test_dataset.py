from unittest import TestCase

from protstruct.dataset import get_pdb_file, extract_chain_data
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


class TestExtractChain(TestCase):
    """
    Test routine for extracting protein chain data
    """

    def test_extract_data(self):
        filename = get_pdb_file( '12as' )
        chain = extract_chain_data( filename, 'A' )

        self.maxDiff = None
        self.assertEqual( len(chain), 328 )

        residue_info = chain[53]
        self.assertEqual(residue_info['aa'], 'V')
        self.assertEqual(residue_info['ss'], 'E')
        self.assertTrue(abs(residue_info['sa'] - 7.0 / 155) <= 1e-2)
        self.assertEqual(residue_info['phi'], -113.0)
        self.assertEqual(residue_info['psi'], 137.0)

        residue_info = chain[299]
        self.assertEqual(residue_info['aa'], 'L')
        self.assertEqual(residue_info['ss'], 'H')
        self.assertTrue(abs(residue_info['sa'] - 2.0 / 170) <= 1e-2)
        self.assertEqual(residue_info['phi'], -76.3)
        self.assertEqual(residue_info['psi'], -38.4)

        tf.gfile.Remove( filename )
