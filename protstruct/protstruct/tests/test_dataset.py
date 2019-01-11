from unittest import TestCase

from protstruct.dataset import get_pdb_file, extract_chain_data, process_pdbs
from protstruct.config import ROOT_DIR, PDB_DIR

import os
import numpy as np
import tensorflow as tf

class TestCreateInputFile(TestCase):
    """
    Test routine for creating input list
    """

    # TODO: test file with wrong pdb id, chain

    def test_create_input_list(self):
        pdb_list = create_input_list(os.path.join(ROOT_DIR, 'data', 'minicullpdb'))
        self.assertEqual(len(pdb_list), 3)

        self.assertEqual(pdb_list[0][0], '12as')
        self.assertEqual(pdb_list[0][1], 'a')

        self.assertEqual(pdb_list[1][0], '16vp')
        self.assertEqual(pdb_list[1][1], 'a')

        self.assertEqual(pdb_list[2][0], '1a0i')
        self.assertEqual(pdb_list[2][1], 'a')

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
        chain = extract_chain_data( filename, 'a' )

        self.maxDiff = None
        self.assertEqual( len(chain), 328 )

        residue_info = chain[53]
        self.assertEqual(residue_info['aa'], 'V')
        self.assertEqual(residue_info['ss'], 'E')
        # see http://prowl.rockefeller.edu/aainfo/volume.htm
        # for residues surface areas
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

