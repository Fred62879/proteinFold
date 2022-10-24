
import csv
import json
import Bio.PDB
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from Bio import SeqIO
from pymol import cmd
from Bio.Seq import Seq
from os.path import join
from Bio.PDB import PDBParser
from biopandas.pdb import PandasPdb
from Bio.SeqRecord import SeqRecord
from Bio.PDB.SASA import ShrakeRupley


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

def extract_residue_from_selection(sel):
    res = []
    #res = []
    #def myfunc(resi,resn,name):
    #    print('%s`%s/%s' % (resn ,resi, name))
    #myspace = {'myfunc': myfunc}
    #cmd.iterate(f'{obj}_interface_R','myfunc(resi,resn,name)', space='myspace') #,res.append((resi,resn))

    objs = cmd.get_object_list(sel)
    for a in range(len(objs)):
        m1 = cmd.get_model(sel + ' and ' + objs[a])
        for x in range(len(m1.atom)):
            if m1.atom[x - 1].resi != m1.atom[x].resi:
                res.append(m1.atom[x].resn)
    return res

def count_num_atoms_for_each_residue(pdb_fn, remove_H=True):
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    if remove_H:
        valid_ids = df['element_symbol'] != 'H'
        df = df[valid_ids]

    num_atoms = []
    chain_ids = read_chain_name_from_pdb(pdb_fn)
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

def read_residue_from_pdb(fn):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    structure = pdb_parser.get_structure("model", fn)

    seq = []
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N','GLY': 'G', 'HIS': 'H','LEU': 'L', 'ARG': 'R', 'TRP': 'W','ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    for chain in structure[0]:
        cur_seq = ''
        for residue in chain:
            res = residue.resname
            if res in d3to1: # ignore hetatm
                cur_seq += d3to1[res]
        seq.append(cur_seq)
    return seq

def read_residue_id_from_pdb(pdb_fn):
    ''' return a list containing unique residue ids for the whole pdb '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    resid_ids = []
    chain_ids = read_chain_name_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        lo, hi = dup_resid_ids[0], dup_resid_ids[-1]
        de_dup_resid_ids = list(np.arange(lo,hi+1))
        resid_ids.extend(de_dup_resid_ids)

    return resid_ids

def read_chain_name_from_pdb(pdb_fn):
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
    chain_ids = read_chain_name_from_pdb(pdb_fn)
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

def find_residue_diff_in_atom_counts(fn1, fn2):
    arr1 = count_num_atoms_for_each_residue(fn1, remove_H=True)
    arr2 = count_num_atoms_for_each_residue(fn2, remove_H=True)
    diff = []
    for map1,map2 in zip(arr1,arr2):
        cur_chain_diff = []
        for k,v1 in map1.items():
            v2 = map2[k]
            if v1 != v2: cur_chain_diff.append([k,v1,v2])
        diff.append(cur_chain_diff)
    return diff

#def prune_pdb_atoms(pdb_fn, out_fn, ranges):
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
    #ppdb.to_pdb(out_fn)
    ppdb.to_pdb(pdb_fn)

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

def set_chain_selector(receptor_chain_id, ligand_chain_id, backbone=False):
    ligand_chain_selector = f'chain {ligand_chain_id}'
    receptor_chain_selector = f'chain {receptor_chain_id}'

    if backbone:
        ligand_chain_selector += ' and backbone'
        receptor_chain_selector += ' and backbone'

    return receptor_chain_selector, ligand_chain_selector

def load_and_select(dist, gt_pdb_fn, pred_pdb_fn, chain_ids, backbone=False, remove_hydrogen=False):
    # Load gt and pred pdb files and select receptor, ligand, and interface
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

    # only deal with dimers
    assert(len(chain_ids) == 2)

    # input processing guarantees first chain is always the receptor chain
    receptor_chain_selector, ligand_chain_selector = set_chain_selector \
        (chain_ids[0], chain_ids[1], backbone)

    for obj in ['native','pred']:
        # select ligand and receptor based on initial assumption
        cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_L', f'{obj} and {ligand_chain_selector}')

        len_r = cmd.count_atoms(f'{obj}_R')
        len_l = cmd.count_atoms(f'{obj}_L')
        assert(len_r > len_l)

        # select receptor interface atoms
        cmd.select(f'{obj}_R_interface_init', f'byres {obj}_R within {dist} of {obj}_L')
        cmd.select(f'interface_{obj}', f'{obj}_R_interface_init + byres {obj}_L within {dist} of {obj}_R_interface_init')

        # color receptor interface in yellow
        cmd.select(f'{obj}_interface_R', f'interface_{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_interface_L', f'interface_{obj} and {ligand_chain_selector}')
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')

def superimpose_receptors(complex_fn):
    # superimpose receptor chains and calculate rmsd for idr, assume existence of corresp sels
    super = cmd.super('native_R','pred_R')
    cmd.color('purple','native_R')
    cmd.color('yellow','native_L')
    cmd.color('gray','pred_R')
    cmd.color('orange','pred_L')
    cmd.multisave(complex_fn, 'all', format='pdb')

def assign_receptor_ligand_chain(pdb_fn, chain_ids):
    ''' Assign two chain ids to be either receptor (long) or ligand (short) '''
    cmd.delete('all')
    cmd.load(pdb_fn, 'pdb')

    # only deal with dimers
    assert(len(chain_ids) == 2)

    # assume 2nd chain is ligand for now
    receptor_chain_selector, ligand_chain_selector = set_chain_selector \
        (chain_ids[0], chain_ids[1])

    # select ligand and receptor based on initial assumption
    cmd.select(f'pdb_R', f'pdb and {receptor_chain_selector}')
    cmd.select(f'pdb_L', f'pdb and {ligand_chain_selector}')

    # make sure ligand has more atoms than receptor, otherwise reverse them
    count_r = cmd.count_atoms(f'pdb_R')
    count_l = cmd.count_atoms(f'pdb_L')
    if count_r < count_l:
        print(f'{pdb_fn[-8:-4]} has receptor and ligand reverted')
        chain_ids = chain_ids[::-1]

def get_sasa(pdb_id, pdb_fn):
    # get solvent accessible surface area of given protein
    p = PDBParser(QUIET=1)
    struct = p.get_structure(pdb_id, pdb_fn)
    sr = ShrakeRupley()
    sr.compute(struct, level="S")
    return struct.sasa

def get_metric_plot_variables(pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn, ranking_fn, chain_id, intrfc_dist):
    ''' Get variables that we will calculate dockq/rmsd value against
        @Return idr_len, num_interface_resid, plddt, sasa, helix chain
    '''

    load_and_select(intrfc_dist, gt_pdb_fn, pred_pdb_fn, chain_id)
    superimpose_receptors(complex_fn)

    # len of ligand
    resid = read_residue_from_pdb(pred_pdb_fn) # 2nd chain is ligand
    len_ligand = len(resid[1])

    # num of interface resid for ligand
    intrfc_resids = extract_residue_from_selection('pred_interface_L')
    len_intrfc_resid = len(intrfc_resids)

    # plddt
    fp = open(ranking_fn)
    rank = json.load(fp)
    plddt = rank["plddts"]["model_1_pred_0"]

    # comprehensive sasa: (sasa_r + sasa_l - sasa_super)/2
    sasa = get_sasa(pdb_id, gt_pdb_fn) + get_sasa(pdb_id, pred_pdb_fn) - get_sasa(pdb_id, complex_fn)
    sasa /= 2

    return (len_ligand, len_intrfc_resid, plddt, sasa)


'''
def load_and_select(dist, gt_pdb_fn, pred_pdb_fn, chain_ids, backbone=False, remove_hydrogen=False):
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

    # only deal with dimers
    assert(len(chain_ids) == 2)

    # assume 2nd chain is ligand for now
    receptor_chain_selector, ligand_chain_selector = set_chain_selector \
        (chain_ids[0], chain_ids[1], backbone)

    for obj in ['native','pred']:
        # select ligand and receptor based on initial assumption
        cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_L', f'{obj} and {ligand_chain_selector}')

        # make sure ligand has more atoms than receptor, otherwise reverse them
        count_r = cmd.count_atoms(f'{obj}_R')
        count_l = cmd.count_atoms(f'{obj}_L')

        if count_r < count_l:
            print(f'{gt_pdb_fn[-8:-4]} has receptor and ligand reverted')
            # update selectors
            receptor_chain_selector, ligand_chain_selector = set_chain_selector \
                (chain_ids[1], chain_ids[0], backbone)

            # redo selection
            cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
            cmd.select(f'{obj}_L', f'{obj} and {ligand_chain_selector}')

        # select receptor interface atoms
        cmd.select(f'{obj}_L_interface_init', f'byres {obj}_L within {dist} of {obj}_L')
        cmd.select(f'interface_{obj}', f'{obj}_L_interface_init + byres {obj}_L within {dist} of {obj}_L_interface_init')

        # color receptor interface in yellow
        cmd.select(f'{obj}_interface_R', f'interface_{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_interface_L', f'interface_{obj} and {ligand_chain_selector}')
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')
'''
