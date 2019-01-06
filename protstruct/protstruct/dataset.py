"""
Utilities to download and process input data to create datasets.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from Bio.PDB import PDBParser, DSSP
from protstruct.utils import download_file
from protstruct.config import PDB_URL, PDB_DIR

import os
import tensorflow as tf

def get_pdb_file( pdb_id, compressed = False ):

    if not tf.gfile.Exists( PDB_DIR ):
        tf.gfile.MakeDirs( PDB_DIR )

    pdb_file = os.path.join( PDB_URL, pdb_id )
    local_pdb_file = os.path.join( PDB_DIR, pdb_id )

    pdb_file += '.pdb'
    local_pdb_file += '.pdb'
    if compressed:
        pdb_file += '.gz'
        local_pdb_file += '.gz'

    # skip if file exists
    if tf.gfile.Exists( local_pdb_file ):
        return local_pdb_file

    return download_file( pdb_file, PDB_DIR )

def extract_chain_data( pdb_file, chain ):
    try:
        end = pdb_file.index('.')
        pdb_id = pdb_file[:end]
    except ValueError:
        pdb_id = pdb_file

    p = PDBParser()
    
    # TODO
    # If using compressed file:
    # UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte
    # - check if compressed, and uncompress if so
    structure = p.get_structure(pdb_id, pdb_file)
    
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp=config.DSSP)

    chain = []
    for dssp_key in filter(lambda t: t[0] == chain, dssp.keys()):
        residue_info = dssp[dssp_key]
        AA, SS, SA, Phi, Psi = residue_info[1:6]
        chain.append( { 'aa':AA, 'ss':SS, 'sa':SA, 'phi':Phi, 'psi':Psi } )

    return chain

#
# TODO
# - logging
# - exception handling, e.g. DSSP may not run for certain PDBs 
# - multiprocessing
#
def process_pdbs(file_list, out):
    with open( game_file, 'r' ) as f:
        for pdb_chain in f:
            # each row is PDB_ID + CHAIN_ID
            pdb_id, chain = pdb_chain[:4], pdb_chain[4]

            pdb_file = get_pdb_file( pdb_id + '.pdb' )
            chain_data = extract_chain_data( pdb_file, chain )

            # TODO: append chain data to file, no with multiprocessing
            
