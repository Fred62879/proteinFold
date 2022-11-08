
import numpy as np
import utils.common as utils

from pymol import cmd
from biopandas.pdb import PandasPdb


def remove_linker(pdb_id, pred_fn, pred_removed_linker_fn, orig_chain_ids, chain_start_ids, gt_chain_bd_ids, n_g, generate_fasta_from_pdb):
    ''' Remove linker from predicted pdb and rename chain
          identifier as per gt pdb
          i) if aa seq comes from gt pdb then also renumber
             residue id as per gt pdb
         ii) if aa seq comes from fasta, then reorder residue id
             such that after removal, residue id is contiguous (1,2,...,n)

        Note: pred pdb residue id is contiguous 1-indexed
          and since all chains are polyg linked, resid id
          is contiguous between chains (polyg linker inclusive)

        Input: orig_chain_ids   identifier for all chain of current gt pdb
               chain_start_ids id of first residue in each chain (from pred pdb)
               gt_chain_bd_ids id of two boundary residue of each chain (from gt pdb)
    '''
    #if exists(pred_removed_linker_fn): return
    n = len(orig_chain_ids)
    assert(n == 2)

    residues = utils.read_residue_from_pdb(pred_fn)[0]
    for i in range(1, len(chain_start_ids)-1):
        id = chain_start_ids[i] - 1 # 1-based to 0-based
        assert(residues[id - n_g : id] == 'G' * n_g)

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pred_fn)
    df = ppdb.df['ATOM']

    acc_linker_length = 0
    linker_atom_ids = []
    atom_los, atom_his, offsets = [], [], []

    for i in range(n): # for each chain
        resid_lo = chain_start_ids[i]
        # resid_hi is residue id of first 'g' in pred
        resid_hi = chain_start_ids[i+1] - n_g

        # make sure gt and pred chain has same length
        # if pred based on fasta, sequence length might differ
        # needs additional processing after linker removal
        gt_resid_lo, gt_resid_hi = gt_chain_bd_ids[i]
        if generate_fasta_from_pdb:
            assert(gt_resid_hi - gt_resid_lo + 1 == resid_hi - resid_lo)

        # locate all atoms in pred pdb for current chain
        atom_ids = df[ (df['residue_number'] >= resid_lo) &
                       (df['residue_number'] < resid_hi)].index
                       #(df['chain_id'] == orig_chain_ids[i]) ].index
        ''' Don't perform chain_id comparison here since linked pdb has
            one single chain. However, this is necessary for gt pdb where
            different chain may have residues with same id.
        '''

        atom_lo, atom_hi = atom_ids[0], atom_ids[-1]

        ''' reorder residue id and rename chain
            df.loc upper bound is inclusive
            atom_hi is id of last atom of current chain
        '''
        if generate_fasta_from_pdb:
            offset = resid_lo - gt_resid_lo
        else:
            offset = acc_linker_length
            acc_linker_length += n_g

        atom_los.append(atom_lo);atom_his.append(atom_hi)
        offsets.append(offset)
        df.loc[atom_lo:atom_hi ,'chain_id'] = orig_chain_ids[i]

        # mark polyg linker for removal
        #if i != n - 2:
        if True:
            linker_atom_lo = atom_hi+1
            linker_resid_hi = chain_start_ids[i+1]
            linker_atom_hi = df[ (df['residue_number'] < linker_resid_hi) ].index[-1]
            linker_atom_ids += list(np.arange(linker_atom_lo, linker_atom_hi+1))

    # re-number chain residue id to be contiguous between chains
    for atom_lo, atom_hi, offset in zip(atom_los, atom_his, offsets):
        df.loc[atom_lo:atom_hi ,'residue_number'] -= offset
    # drop linker residues
    df.drop(linker_atom_ids, axis=0, inplace=True)

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(pred_removed_linker_fn)

def remove_extra_residue_predicted_pdb(in_fn, out_fn, gt_pdb_fn, orig_chain_ids, gt_chain_bd_ids):
    seqs1 = utils.read_residue_from_pdb(in_fn)
    seqs2 = utils.read_residue_from_pdb(gt_pdb_fn)
    ranges = find_prune_ranges_all_chains(seqs1, seqs2, orig_chain_ids)
    prune_seq_given_ranges(in_fn, out_fn, orig_chain_ids, ranges[0], gt_chain_bd_ids)

def renumber_residue_perdicted_pdb(in_fn, out_fn, orig_chain_ids, gt_chain_bd_ids):
    ''' renumber predicted pdb residue to be identical to gt pdb ids '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fn)
    df = ppdb.df['ATOM']

    renumber_atom_ids, renumber_offsets = [], []
    for chain_id, gt_chain_bd_id in zip(orig_chain_ids, gt_chain_bd_ids):
        ids = df.loc[ (df['chain_id'] == chain_id) ].index
        resid_lo = df.loc[ids[0],'residue_number']
        offset = resid_lo - gt_chain_bd_id[0]
        renumber_atom_ids.append(ids)
        renumber_offsets.append(offset)

    for ids, renumber_offset in zip(renumber_atom_ids, renumber_offsets):
        df.loc[ids, 'residue_number'] -= renumber_offset

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def calculate_rmsd(pdb_id, gt_pdb_fn, pred_pdb_fn, interface_dist, order_chain_ids, remove_backbone, remove_hydrogen, verbose=False):
    # calculate rmsd between gt and pred via superimposing pred onto gt (with receptor only)
    rmsds = []
    utils.load_and_select \
        (interface_dist, gt_pdb_fn, pred_pdb_fn, ordered_chain_ids,
         backbone=backbone, remove_hydrogen=remove_hydrogen)

    # superimpose receptor chains and calculate rmsd for ligand
    cmd.super('native_R','pred_R')
    rmsd = cmd.rms_cur('native_L','pred_L')
    rmsds.append(rmsd)

    # save two objects after superimposing receptor chain
    cmd.color('gray','native')
    cmd.color('red','pred')
    for obj in ['native','pred']:
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')

    # calculate rmsd for interface idr only
    rmsd = cmd.rms_cur('native_interface_L','pred_interface_L')
    rmsds.append(rmsd)
    rmsds = [round(rmsd,3) for rmsd in rmsds]
    return rmsds

################
# residue prune
def find_subseq_range(seq1, seq2):
    n, m = len(seq1), len(seq2)

    # assume seq1 has extra residues at head
    # e.g. seq1: ABCD, seq2: BCD...
    for i in range(n):
        if m > n - i:
            # seq2 has extra residues at tail
            # e.g. seq1: ABCD, seq2: BCDE
            if seq1[i:] == seq2[:n-i]:
                return np.array([[0,i]]), np.array([[n-i,m]])
        else: # seq1 has extra residues at tail
            if seq1[i:i+m] == seq2:
                return np.array([[0,i],[i+m,n]]), []

    # assume seq2 has extra residues at head
    for i in range(m):
        if n > m - i:
            # seq1 has extra residues at tail
            if seq1[:m-i] == seq2[i:]:
                return np.array([[m-i,n]]), np.array([[0,i]])
        else: # seq2 has extra residues at tail
            if seq2[i:i+n] == seq1:
                return [], np.array([[0,i],[i+n,m]])
    raise Exception('cannot find subseq of pred and native that match')

def find_prune_ranges_all_chains(seqs1, seqs2, chain_ids):
    ''' Find range of residues to prune for each chain of the two sequences
        which will be identical after pruning.
        Assume the two sequences differ only in residues at two ends.
          e.g. seqs1 'BCDE' seqs2 'ABCD' where core residues 'BC' are
               shared and only residues at the two ends differ.
        Returned range is residue id (1-based)
    '''
    rnge1, rnge2 = [], []
    acc_len1, acc_len2 = 0, 0
    for seq1, seq2, chain_id in zip(seqs1, seqs2, chain_ids):
        cur_rnge1, cur_rnge2 = find_subseq_range(seq1, seq2)
        if len(cur_rnge1) != 0:
            rnge1.append(cur_rnge1 + acc_len1)
        if len(cur_rnge2) != 0:
            rnge2.append(cur_rnge2 + acc_len2)
        acc_len1 += len(seq1)
        acc_len2 += len(seq2)
    return np.array(rnge1)+1, np.array(rnge2)+1

def prune_seq_given_range(df, rnge, chain_id, to_remove_atom_ids, renumber_atom_ids, renumber_offsets, prev_hi, acc_offset):
    ''' Remove atoms falling within residue id lo (inclusive) and hi
          (exclusive) and renumber to make sure after pruning, residue
          id goes from 1-n consecutively
    '''

    (resid_lo,resid_hi) = rnge
    atom_ids = df[ (df['chain_id'] == chain_id) &
                   (df['residue_number'] >= resid_lo) &
                   (df['residue_number'] < resid_hi) ].index
    to_remove_atom_ids.extend(atom_ids)

    if prev_hi != -1:
        renumber_range = [prev_hi, resid_lo]
        renumber_offsets.append(acc_offset)

        # update residue number (1-based)
        cur_atom_ids = df[ (df['chain_id'] == chain_id) &
                           (df['residue_number'] >= prev_hi) &
                           (df['residue_number'] < resid_lo) ].index
        renumber_atom_ids.append(cur_atom_ids)

    acc_offset += (resid_hi - resid_lo)
    return resid_hi, acc_offset

def prune_seq_given_ranges(in_fn, out_fn, chain_ids, ranges, gt_chain_bd_ids):
    ''' Prune extra residues from predicted pdb so that it match with gt
        Also renumber residue id to be the same as in gt
        Only used if aa seq comes from fasta file
    '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fn)
    df = ppdb.df['ATOM']

    prev_hi, acc_offset = -1, 0
    to_remove_atom_ids = []
    renumber_atom_ids, renumber_offsets = [], []

    # gather atom ids for extra residues &
    # atom ids for renumbering with corresponding offset
    for i, (cur_ranges, chain_id) in enumerate(zip(ranges, chain_ids)):
        for rnge in cur_ranges:
            if rnge[0] == rnge[1]: continue
            prev_hi, acc_offset = prune_seq_given_range\
                (df, rnge, chain_id, to_remove_atom_ids, renumber_atom_ids,
                 renumber_offsets, prev_hi, acc_offset)

    # remove extra residues
    df.drop(to_remove_atom_ids, axis=0, inplace=True)
    # renumber such that residue id is contiguous (1,2,...n)
    for ids, renumber_offset in zip(renumber_atom_ids, renumber_offsets):
        df.loc[ids, 'residue_number'] -= renumber_offset

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)
