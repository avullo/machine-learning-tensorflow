"""
Utilities to download and process input data to create datasets.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from Bio.PDB import PDBParser, DSSP
from protstruct.utils import download_file
import protstruct.config as config # import PDB_URL, PDB_DIR, DSSP

import os
import h5py
import tensorflow as tf

import pprint

def create_input_list(pdb_list_fname):
    """
    create a list of tuples (pdb_id, chain) from a text file
    """
    pdb_list = []
    with open(pdb_list_fname, 'r') as f:
        for record in f.read().splitlines():
            pdb_id, chain = record[:-1], record[-1]

            # check PDB ID and chain are valid
            if not pdb_id.isalnum() or len(pdb_id) != 4 or not chain.isalpha() or len(chain) != 1:
                continue

            pdb_list.append((pdb_id, chain))

    return pdb_list

def get_pdb_file( pdb_id, compressed = False ):
    pdb_id = pdb_id.lower()
    
    if not tf.gfile.Exists( config.PDB_DIR ):
        tf.gfile.MakeDirs( config.PDB_DIR )

    pdb_file = os.path.join( config.PDB_URL, pdb_id )
    local_pdb_file = os.path.join( config.PDB_DIR, pdb_id )

    pdb_file += '.pdb'
    local_pdb_file += '.pdb'
    if compressed:
        pdb_file += '.gz'
        local_pdb_file += '.gz'

    # skip if file exists
    if tf.gfile.Exists( local_pdb_file ):
        return local_pdb_file

    return download_file( pdb_file, config.PDB_DIR )

def extract_chain_data( pdb_file, chain ):
    try:
        end = os.path.basename(pdb_file).index('.')
        pdb_id = os.path.basename(pdb_file)[:end]
    except ValueError:
        pdb_id = os.path.basename(pdb_file)

    p = PDBParser()
    
    # TODO
    # If using compressed file:
    # UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte
    # - check if compressed, and uncompress if so
    structure = p.get_structure(pdb_id, pdb_file)
    
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    chain_data = []
    chain = chain.upper() # DSSP records chains in upper case (TODO: always?)

    for dssp_key in filter(lambda t: t[0] == chain, dssp.keys()):
        residue_info = dssp[dssp_key]
        AA, SS, SA, Phi, Psi = residue_info[1:6]
        chain_data.append( { 'AA':AA, 'SS':SS, 'SA':SA, 'Phi':Phi, 'Psi':Psi } )

    return chain_data

#
# TODO
# - logging
# - exception handling:
#   - DSSP may not run for certain PDBs
#   - discontinuous chain (PDBConstructionWarning)
# - multiprocessing?
#
def process_pdbs(pdb_list, out_filename, min_chain_len=30):
    with h5py.File(out_filename, 'w') as f:
        for pdb_id, chain in pdb_list:
            try:
                # chain data dedcated group
                grp = f.create_group(pdb_id + '/' + chain)
                data = extract_chain_data( get_pdb_file( pdb_id ), chain )

                # filter out short chains
                if len(data) < min_chain_len:
                    grp.attrs['status'] = 0
                    continue

                # special type to write variable length strings to HDF5
                dt = h5py.special_dtype(vlen=bytes)

                # AA/SS/SA/Phi/Psi datasets
                # TODO: phi/psi datasets
                grp.create_dataset('AA', (1,), dtype=dt, data=''.join(map(lambda r: r['AA'], data)))
                grp.create_dataset('SS', (1,), dtype=dt, data=''.join(map(lambda r: r['SS'], data)))
                grp.create_dataset('SA',  data=list(map(lambda r: r['SA'], data)))
                grp.create_dataset('Phi', data=list(map(lambda r: r['Phi'], data)))
                grp.create_dataset('Psi', data=list(map(lambda r: r['Psi'], data)))

                grp.attrs['status'] = 1

            except:
                grp.attrs['status'] = 0
            
        # TODO: append chain data to file, no with multiprocessing
