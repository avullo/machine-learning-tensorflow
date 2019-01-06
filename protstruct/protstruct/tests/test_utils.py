import os
import tensorflow as tf

from unittest import TestCase
from protstruct.utils import download_file

from protstruct.config import PDB_DIR

class TestDownloadFile(TestCase):
    """
    Test routine for downloading files
    """

    def test_download_file_no_dest(self):
        filename = download_file( 'https://files.rcsb.org/download/12as.pdb.gz' )
        
        self.assertEqual( filename, '12as.pdb.gz' )
        self.assertTrue( tf.gfile.Exists( filename ) )

        tf.gfile.Remove( filename )

    def test_download_file_dest(self):
        filename = download_file( 'https://files.rcsb.org/download/12as.pdb.gz', PDB_DIR )
        
        self.assertEqual( filename, os.path.join( PDB_DIR, '12as.pdb.gz' ) )
        self.assertTrue( tf.gfile.Exists( filename ) )

        tf.gfile.Remove( filename )

