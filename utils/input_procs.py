
import pickle
import numpy as np
import utils.common as utils

from os import listdir
from os.path import exists, join
from biopandas.pdb import PandasPdb


def poly_g_link_fasta(indir, outdir, chain_start_ids, fasta_group, poly_g, n_g):
    ''' poly glycine link all chains of a given protein where sequence
          of each chain is downloaded from rcsb.org
        @Param
          fasta_group: name of fasta files, each containing one chain of the cur protein
    '''
    multimer = ''
    cur_chain_start_ids = []
    pdb_id = (fasta_group[0]).split('_')[0]

    for fn in fasta_group:
        monomer = utils.read_fasta(join(indir, fn + '.fasta'))
        cur_chain_start_ids.append(len(multimer)+1) # 1-based
        multimer += monomer + poly_g

    # add last start_id for convenience of polyg removal (for loop)
    cur_chain_start_ids.append(len(multimer) + 1) # resid id 1-based
    multimer = multimer[:-n_g]

    ofn = join(outdir, pdb_id + '.fasta')
    utils.save_fasta(multimer, ofn, id=pdb_id)
    chain_start_ids[pdb_id] = cur_chain_start_ids

def poly_g_link_pdb(seqs, n_g, input_fasta_dir, chain_start_ids):
    ''' poly glycine link all chains of a given protein where sequence
          of each chain is obtained from the corspd. pdb file
    '''
    fasta = [reduce(lambda acc, seq: acc + seq + linker, seqs, '')[:-n_g]]
    out_fn = join(input_fasta_dir, id + '.fasta')
    utils.save_fasta(fasta, out_fn)

    # get id of start residues for each chain
    acc, start_ids = 1, []
    for seq in seqs:
        start_ids.append(acc)
        acc += len(seq) + n_g
    start_ids.append(len(fasta) + n_g + 1) # for ease of removal
    chain_start_ids[id] = start_ids

def combine_source_fasta(source_fasta_dir):
    # combine fasta files (each contain one chain) to a single fasta file
    prev_id = ''
    cur_group, groups = [], []

    fns = sorted(listdir(source_fasta_dir))
    for fn in fns:
        if '.fasta' not in fn: continue
        cur_id = fn.split('_')[0]
        if cur_id != prev_id:
           groups.append(cur_group)
           cur_group = [fn.split('.')[0]]
           prev_id = cur_id
        else:
           cur_group.append(fn.split('.')[0])
    groups.append(cur_group)
    return groups[1:]

def poly_g_link_all(strategy, source_fasta_dir, input_fasta_dir, chain_start_resid_ids_fn, n_g):
    if exists(chain_start_resid_ids_fn) and exists(input_fasta_dir):
        pdb_ids1 = utils.parse_pdb_ids(source_fasta_dir, '.fasta')
        pdb_ids2 = utils.parse_pdb_ids(input_fasta_dir, '.fasta')
        bypass = set(pdb_ids1) == set(pdb_ids2)
    else: bypass = False
    if bypass: return

    fasta_groups = combine_source_fasta(source_fasta_dir)
    print('= poly g linking fasta')
    poly_g = 'G' * n_g
    chain_start_resid_ids = {}

    for fasta_group in fasta_groups:
        poly_g_link_fasta(source_fasta_dir, input_fasta_dir,
                          chain_start_resid_ids, fasta_group, poly_g, n_g)
    with open(chain_start_resid_ids_fn, 'wb') as fp:
        pickle.dump(chain_start_resid_ids, fp)

def prune_pdb_atoms(pdb_fn, ranges):
    ''' Remove extra atoms from pdb files and renumber atom ids
        @Param
          ranges (upper exclusive)
    '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    # offset between atom id and dataframe id
    # increment by 1 after each chain
    atom_offset = 0

    id_lo, acc_len = 0, 0
    to_remove_atom_ids = []
    acc_lens, decrmt_rnge = [], []

    # gather ranges of atoms to delete and to decrement with id
    for i, cur_chain_ranges in enumerate(ranges):
        for rnge in cur_chain_ranges:
            (lo, hi) = rnge
            ids = list(np.arange(lo-atom_offset, hi+1-atom_offset))
            to_remove_atom_ids.extend(ids)

            # update range of atoms that are kept
            decrmt_rnge.append([id_lo, lo - 1]) #-acc_len
            id_lo = lo #-acc_len

            # update acc_len
            acc_lens.append(acc_len)
            acc_len += hi - lo + 1
        atom_offset += 1

    acc_lens.append(acc_len)
    decrmt_rnge.append([id_lo, len(df)-1])

    # drop atoms
    df.drop(to_remove_atom_ids, axis=0, inplace=True)

    # decrement atoms
    for i, (length, rnge) in enumerate(zip(acc_lens, decrmt_rnge)):
        lo, hi = rnge
        #print(lo, hi, length)
        df.loc[lo:hi,'atom_number'] -= length # upper inclusive

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(pdb_fn)
