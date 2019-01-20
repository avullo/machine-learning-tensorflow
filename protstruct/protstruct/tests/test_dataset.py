from unittest import TestCase

from protstruct.dataset import * # get_pdb_file, extract_chain_data, process_pdbs
from protstruct.config import ROOT_DIR, PDB_DIR

import os
import tempfile
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

        self.assertEqual(pdb_list[0][0], '12AS')
        self.assertEqual(pdb_list[0][1], 'A')

        self.assertEqual(pdb_list[1][0], '16VP')
        self.assertEqual(pdb_list[1][1], 'A')

        self.assertEqual(pdb_list[2][0], '1A0I')
        self.assertEqual(pdb_list[2][1], 'A')

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
        self.assertEqual(residue_info['AA'], 'V')
        self.assertEqual(residue_info['SS'], 'E')
        # see http://prowl.rockefeller.edu/aainfo/volume.htm
        # for residues surface areas
        self.assertTrue(abs(residue_info['SA'] - 7.0 / 155) <= 1e-2)
        self.assertEqual(residue_info['Phi'], -113.0)
        self.assertEqual(residue_info['Psi'], 137.0)

        residue_info = chain[299]
        self.assertEqual(residue_info['AA'], 'L')
        self.assertEqual(residue_info['SS'], 'H')
        self.assertTrue(abs(residue_info['SA'] - 2.0 / 170) <= 1e-2)
        self.assertEqual(residue_info['Phi'], -76.3)
        self.assertEqual(residue_info['Psi'], -38.4)

        tf.gfile.Remove( filename )

class TestProcessPDBs(TestCase):
    """
    Test routine to process multiple PDB chains
    """

    def test_process_pdbs(self):
        pdb_list = create_input_list(os.path.join(ROOT_DIR, 'data', 'minicullpdb'))

        h5fname = os.path.join(ROOT_DIR, 'data', next(tempfile._get_candidate_names()) + '.h5')
        process_pdbs(pdb_list, h5fname)
        self.assertTrue( tf.gfile.Exists( h5fname ) )

        # check file content
        with h5py.File(h5fname, 'r') as f:
            with self.assertRaises(Exception) as context:
                f['12AS/B']

            chain = f['12AS/A']
            aa = np.array(chain['AA'])[0].decode('ascii')
            ss = np.array(chain['SS'])[0].decode('ascii')
            sa = np.array(chain['SA'])

            self.assertEqual(aa, 'AYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLLN')
            self.assertEqual(ss, '-HHHHHHHHHHHHHHHHHHHHHHH-EEE----SEEETTSS-S--TTTT----EE--SSSTT--EEE-S--TTHHHHHHHHTT--TT-EEEEEEEEE-TT-S---SS--SEEEEEEEEEE--TT--SHHHHHHHHHHHHHHHHHHHHHHHHHS-----S-SS-EEEEHHHHHHHSSSS-HHHHHHHHHHHHSEEEEE--SS--SSS--SS---TTTB--SSB-TTSSB-SEEEEEEEETTTTEEEEEEEEEEB--HHHHHHHHHHHT-TTGGGSHHHHHHHTT-S--EEEEEEEHHHHHHHHHT-S-GGGTS-----HHHHHH--S---')

            self.assertEqual(len(sa), 328)
            try:
                np.testing.assert_allclose(sa[9:19], [0.23076923076923078,0.27918781725888325,0.0,0.13658536585365855,0.6153846153846154,0.16847826086956522,0.0,0.07692307692307693,0.5241935483870968,0.08080808080808081], rtol=1e-2)
                np.testing.assert_allclose(sa[-2:], [0.35365853658536583,0.09554140127388536], rtol=1e-2)
            except AssertionError:
                self.fail('Failed SA for 12AS/A')

            # TODO: test Phi/Psi values

        tf.gfile.Remove( h5fname )
