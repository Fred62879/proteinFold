
import csv
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from Bio import SeqIO
from pymol import cmd
from Bio.Seq import Seq
from os.path import join
from biopandas.pdb import PandasPdb
from Bio.SeqRecord import SeqRecord

def read_line_by_line(fn):
    f = open(fn, 'r')
    lines = f.readlines()
    return [line.strip() for line in lines]

def write_to_csv(data, out_fn):
    with open(out_fn, 'w', newline='') as fp:
        writer = csv.writer(fp)
        for row in data:
            writer.writerow(row)
    #np.savetxt(args.rmsd_fn, np.array(rmsds),delimiter=',',
    #           header=args.rmsd_names,comments='')

def read_fasta(fn):
    fasta_sequences = SeqIO.parse(open(fn),'fasta')
    for fasta in fasta_sequences:
        monomer = str(fasta.seq)
        break
    return monomer

def save_fasta(seq, fn, id='0'):
    seq = SeqRecord(Seq(seq),id=id)
    SeqIO.write(seq, fn, 'fasta')

def read_num_atoms_for_each_residue(pdb_fn, remove_H=True):
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    if remove_H:
        valid_ids = df['element_symbol'] != 'H'
        df = df[valid_ids]

    num_atoms = []
    chain_ids = read_chain_id_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_num_atoms = {}
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        lo, hi = dup_resid_ids[0], dup_resid_ids[-1]

        for resid_id in range(lo, hi+1):
            valid_atoms = df[ (df['chain_id']== chain_id) &
                              (df['residue_number'] == resid_id) ]
            n_atoms = len(valid_atoms)
            cur_num_atoms[resid_id] = n_atoms
        num_atoms.append(cur_num_atoms)
    return num_atoms

def read_residue_from_pdb(fn, print=False):
    '''
    seq = ''
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N','GLY': 'G', 'HIS': 'H','LEU': 'L', 'ARG': 'R', 'TRP': 'W','ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    for model in structure:
        for chain in model:
            for residue in chain:
                res = residue.resname
                if res in d3to1:
                    seq += d3to1[res]
        return seq
    '''
    fastas = []
    with open(fn, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            if print:
                print('>' + record.id)
                print(record.seq)
            else:
                _ = record.id
                fastas.append(str(record.seq))
    return fastas

def read_residue_id_from_pdb(pdb_fn):
    ''' return a list containing unique residue ids for the whole pdb '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    resid_ids = []
    chain_ids = read_chain_id_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        lo, hi = dup_resid_ids[0], dup_resid_ids[-1]
        de_dup_resid_ids = list(np.arange(lo,hi+1))
        resid_ids.extend(de_dup_resid_ids)

    return resid_ids

def read_chain_id_from_pdb(pdb_fn):
    ''' read chain identifier, de-depulicate while preserving order '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']
    return list(dict.fromkeys(df.loc[:,'chain_id']))

def read_chain_bd_resid_id_from_pdb(pdb_fn, pdb_id, bd_resid_ids):
    ''' Find residue id of two boundary residues in each chain '''
    print(f'Reading chain boundary residue id for {pdb_id}')

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']
    dfh = ppdb.df['HETATM']

    cur_bd_resid_ids = []
    chain_ids = read_chain_id_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        resid_lo, resid_hi = dup_resid_ids[0], dup_resid_ids[-1]

        '''
        # to keep resid id for hetatm
        cur_chain_htm_atom_ids =dfh['chain_id'] == chain_id
        if len(cur_chain_htm_atom_ids) != 0:
            dup_resid_ids = list(dfh.loc[cur_chain_htm_atom_ids,'residue_number'])
            if len(dup_resid_ids) != 0:
                htm_resid_lo, htm_resid_hi = dup_resid_ids[0], dup_resid_ids[-1]
                print(resid_lo, resid_hi, htm_resid_lo, htm_resid_hi)
                if resid_lo > htm_resid_lo:
                    print('HETATOM lower')
                if resid_hi < htm_resid_hi:
                    print('HETATOM higher')
                resid_lo = min(resid_lo, htm_resid_lo)
                resid_hi = max(resid_hi, htm_resid_hi)
        '''

        cur_bd_resid_ids.append([resid_lo, resid_hi])

    bd_resid_ids[pdb_id] = cur_bd_resid_ids
    return bd_resid_ids

def poly_g_link(indir, outdir, chain_start_ids, fasta_group, poly_g, n_g):
    multimer = ''
    cur_chain_start_ids = []
    pdb_id = (fasta_group[0]).split('_')[0]

    for fn in fasta_group:
        monomer = read_fasta(join(indir, fn + '.fasta'))
        cur_chain_start_ids.append(len(multimer)+1) # 1-based
        multimer += monomer + poly_g

    # add last start_id for convenience of polyg removal (for loop)
    cur_chain_start_ids.append(len(multimer) + 1) # resid id 1-based
    multimer = multimer[:-n_g]

    ofn = join(outdir, pdb_id + '.fasta')
    save_fasta(multimer, ofn, id=pdb_id)
    chain_start_ids[pdb_id] = cur_chain_start_ids

def find_subseq_range(seq1, seq2):
    n, m = len(seq1), len(seq2)
    for i in range(n):
        if seq1[i:i+m] == seq2:
            return np.array([[0,i],[i+m,n]])
    assert(False)

def find_residue_diff_in_atom_counts(fn1, fn2):
    arr1 = read_num_atoms_for_each_residue(fn1, remove_H=True)
    arr2 = read_num_atoms_for_each_residue(fn2, remove_H=True)
    diff = []
    for map1,map2 in zip(arr1,arr2):
        cur_chain_diff = []
        for k,v1 in map1.items():
            v2 = map2[k]
            if v1 != v2: cur_chain_diff.append([k,v1,v2])
        diff.append(cur_chain_diff)
    return diff

def prune_extra_atoms(pdb_fn, out_fn, ranges):
    ''' remove extra atoms from pdb files and renumber atom ids '''
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
    ppdb.to_pdb(out_fn)

def find_prune_ranges_all_chains(seqs1, seqs2, chain_ids, prune_X=False):
    ''' assume seq in seqs1 is longer than corresponding seq in seq2
        returned range is residue id (1-based)
    '''
    acc_len = 0
    rnge1, rnge2 = [], []
    for seq1, seq2, chain_id in zip(seqs1, seqs2, chain_ids):
        if prune_X:
            cur_rnge1 = find_x_range(seq1)
            if len(cur_rnge1) != 0:
                rnge1.append(cur_rnge1 + acc_len)
            cur_rnge2 = find_x_range(seq2)
            if len(cur_rnge2) != 0:
                rnge2.append(cur_rnge2 + acc_len)
        else:
            cur_rnge = find_subseq_range(seq1, seq2)
            if len(cur_rnge) != 0:
                rnge1.append(cur_rnge + acc_len)
        acc_len += len(seq1)

    return np.array(rnge1)+1, np.array(rnge2)+1

def prune_renumber_seq_given_range(df, rnge, chain_id, to_remove_atom_ids, renumber_atom_ids, renumber_offsets, prev_hi, acc_offset):
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

def renumber_seq_per_gt(df, chain_ids, gt_chain_bd_ids, renumber_atom_ids, renumber_offsets):
    for chain_id, gt_chain_bd_id in zip(chain_ids, gt_chain_bd_ids):
        ids = df.loc[ (df['chain_id'] == chain_id) ].index
        resid_lo = df.loc[ids[0],'residue_number']
        offset = resid_lo - gt_chain_bd_id[0]
        renumber_atom_ids.append(ids)
        renumber_offsets.append(offset)

def prune_renumber_seq_given_ranges(in_fn, out_fn, chain_ids, ranges, gt_chain_bd_ids):
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
            if rnge[0] == rnge[1]:
                continue
            prev_hi, acc_offset = prune_renumber_seq_given_range\
                (df, rnge, chain_id, to_remove_atom_ids, renumber_atom_ids,
                 renumber_offsets, prev_hi, acc_offset)

    # remove extra residues
    df.drop(to_remove_atom_ids, axis=0, inplace=True)

    # renumber such that residue id is contiguous (1,2,...n)
    for ids, renumber_offset in zip(renumber_atom_ids, renumber_offsets):
        df.loc[ids, 'residue_number'] -= renumber_offset

    # renumber such that residue id match gt
    renumber_atom_ids, renumber_offsets = [], []
    renumber_seq_per_gt(df, chain_ids, gt_chain_bd_ids, renumber_atom_ids, renumber_offsets)
    for ids, renumber_offset in zip(renumber_atom_ids, renumber_offsets):
        df.loc[ids, 'residue_number'] -= renumber_offset

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def load_and_select(dist, gt_pdb_fn, pred_pdb_fn, chain_names, backbone=False, remove_hydrogen=False):
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

    idr_chain_id = chain_names[-1]
    target_chain_selector = f'chain {idr_chain_id}'
    receptor_chain_selector = f'not chain {idr_chain_id}'

    if backbone:
        target_chain_selector += ' and backbone'
        receptor_chain_selector += ' and backbone'

    for obj in ['native','pred']:
        cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_T', f'{obj} and {target_chain_selector}')
        cmd.select(f'{obj}_R_interface_init', f'byres {obj}_R within {dist} of {obj}_T')
        cmd.select(f'interface_{obj}', f'{obj}_R_interface_init + byres {obj}_T within {dist} of {obj}_R_interface_init')

        # color receptor interface of {obj} in yellow
        cmd.select(f'{obj}_interface_R', f'interface_{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_interface_T', f'interface_{obj} and {target_chain_selector}')
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_T')
