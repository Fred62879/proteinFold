
import csv
import json
import subprocess
import numpy as np

import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

from Bio import SeqIO
from pymol import cmd
from os import listdir
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from os.path import join, exists
from biopandas.pdb import PandasPdb
from Bio.SeqRecord import SeqRecord
from Bio.PDB.SASA import ShrakeRupley


def run_bash(cmd):
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def read_csv(fn):
    f = open(fn, 'r')
    lines = f.readlines()
    return [line.strip() for line in lines]

def parse_csv(fn):
    pdb_ids, res = [], {}
    data = read_csv(fn)
    metric_names = data[0].split(',')[1:]
    for i in range(1, len(data)):
        entries = data[i].split(',')
        pdb_id = entries[0]
        pdb_ids.append(pdb_id)
        cur_data = {}
        for i, name in enumerate(metric_names):
            cur_data[name] = float(entries[i + 1])
        res[pdb_id] = cur_data
    return pdb_ids, metric_names, res

def plot_scatter(dir, data, pdb_ids):
    stat_names = ['idr_len','num_interface_resid','plddt','sasa']
    metric = []
    for pdb_id in pdb_ids:
        metric.append(data[pdb_id]['dockQ'])

    for stat_name in stat_names:
        cur_stat = []
        for pdb_id in pdb_ids:
            cur_stat.append(data[pdb_id][stat_name])

        plt.scatter(cur_stat, metric)
        plt.xlabel(stat_name)
        plt.ylabel('dockQ')
        plt.savefig(join(dir, stat_name + '.png'))
        plt.close()

def write_to_csv(data, out_fn):
    with open(out_fn, 'w', newline='') as fp:
        writer = csv.writer(fp)
        for row in data:
            writer.writerow(row)
    #np.savetxt(args.rmsd_fn, np.array(rmsds),delimiter=',',
    #           header=args.rmsd_names,comments='')

################
# residue procs
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

def read_residue_from_pdb(fn):
    pdb_parser = PDBParser(QUIET = True)
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
    #print(f'Reading chain boundary residue id for {pdb_id}')

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

def assert_seq_equality(gt_pdb_fn, pred_pdb_fn):
    # check if sequence from processed prediction pdb match sequence from gt pdb
    gt_pdb = read_residue_from_pdb(gt_pdb_fn)
    pred_pdb = read_residue_from_pdb(pred_pdb_fn)
    if not gt_pdb == pred_pdb:
        pdb_id = gt_pdb_fn[:4]
        raise Exception(f'{pdb_id}: prediction and gt pdb sequence don\'t match')

def assert_equal_num_chains(gt_pdb_fn, pred_pdb_fn):
    gt_pdb = read_residue_from_pdb(gt_pdb_fn)
    pred_pdb = read_residue_from_pdb(pred_pdb_fn)
    if not len(gt_pdb) == len(pred_pdb):
        pdb_id = gt_pdb_fn[:4]
        raise Exception(f'{pdb_id}: prediction and gt pdb have diff num of chains')

##############
# fasta procs
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
        monomer = read_fasta(join(indir, fn + '.fasta'))
        cur_chain_start_ids.append(len(multimer)+1) # 1-based
        multimer += monomer + poly_g

    # add last start_id for convenience of polyg removal (for loop)
    cur_chain_start_ids.append(len(multimer) + 1) # resid id 1-based
    multimer = multimer[:-n_g]

    ofn = join(outdir, pdb_id + '.fasta')
    save_fasta(multimer, ofn, id=pdb_id)
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
        pdb_ids1 = parse_pdb_ids(source_fasta_dir, '.fasta')
        pdb_ids2 = parse_pdb_ids(input_fasta_dir, '.fasta')
        bypass = set(pdb_ids1) == set(pdb_ids2)
    else: bypass = False
    if bypass: return

    fasta_groups = combine_source_fasta(source_fasta_dir)
    print('= poly g linking fasta')
    poly_g = 'G' * n_g
    chain_start_resid_ids = {}

    for fasta_group in fasta_groups:
        utils.poly_g_link(source_fasta_dir, input_fasta_dir,
                          chain_start_resid_ids, fasta_group, poly_g, n_g)
    with open(self.chain_start_resid_ids_fn, 'wb') as fp:
        pickle.dump(chain_start_resid_ids, fp)

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

#####################
# prune and renumber
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


#########
# others
def parse_pdb_ids(dir, suffix):
    # collect all unique pdb ids in a given directory
    #pdb_ids = next(os.walk(dir))[1]
    pdb_ids = list(set(listdir(dir)))
    if len(pdb_ids) == 0: return []
    pdb_ids = [pdb_id[:4] for pdb_id in pdb_ids if suffix in pdb_id]
    return np.array(list(set(pdb_ids)))

def set_chain_selector(receptor_chain_id, ligand_chain_id, backbone=False):
    ligand_chain_selector = f'chain {ligand_chain_id}'
    receptor_chain_selector = f'chain {receptor_chain_id}'

    if backbone:
        ligand_chain_selector += ' and backbone'
        receptor_chain_selector += ' and backbone'

    return receptor_chain_selector, ligand_chain_selector

def load_and_select(interface_dist, gt_pdb_fn, pred_pdb_fn, ordered_chain_ids, backbone=False, remove_hydrogen=False):
    # Load gt and pred pdb files and select receptor, ligand, and interface
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

    # only deal with dimers
    assert(len(ordered_chain_ids) == 2)

    # first chain is receptor
    receptor_chain_selector, ligand_chain_selector = set_chain_selector \
        (ordered_chain_ids[0], ordered_chain_ids[1], backbone)

    for obj in ['native','pred']:
        # select ligand and receptor based on initial assumption
        cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_L', f'{obj} and {ligand_chain_selector}')

        # select receptor interface atoms
        cmd.select(f'{obj}_R_interface_init', f'byres {obj}_R within {interface_dist} of {obj}_L')
        cmd.select(f'interface_{obj}', f'{obj}_R_interface_init + byres {obj}_L within {interface_dist} of {obj}_R_interface_init')

        # color receptor interface in yellow
        cmd.select(f'{obj}_interface_R', f'interface_{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_interface_L', f'interface_{obj} and {ligand_chain_selector}')
        cmd.color('yellow', f'{obj}_interface_R')
        cmd.color('blue',f'{obj}_interface_L')

def assign_receptor_ligand(pdb_fn, chain_ids):
    cmd.delete('all')
    cmd.load(pdb_fn, 'native')

    # assume first chain is receptor
    (r, l) = chain_ids
    reverted = False
    receptor_chain_selector, ligand_chain_selector = set_chain_selector(r, l)

    # select ligand and receptor based on initial assumption
    cmd.select(f'native_R', f'native and {receptor_chain_selector}')
    cmd.select(f'native_L', f'native and {ligand_chain_selector}')

    # make sure ligand has less atoms than receptor, otherwise revert them
    count_r = cmd.count_atoms(f'native_R')
    count_l = cmd.count_atoms(f'native_L')
    if count_r < count_l:
        #print(f'{pdb_fn[-8:-4]} has receptor and ligand reverted')
        (l, r) = chain_ids
        reverted = True

    return (r, l, reverted)

def superimpose_receptors():
    # superimpose receptor chains and calculate rmsd for idr, assume existence of corresp sels
    super = cmd.super('native_R','pred_R')
    #cmd.color('purple','native_R')
    #cmd.color('yellow','native_L')
    #cmd.color('gray','pred_R')
    #cmd.color('orange','pred_L')
    #cmd.multisave(complex_fn, 'all', format='pdb')

def get_sasa(pdb_id, pdb_fn):
    # get solvent accessible surface area of given protein
    p = PDBParser(QUIET=1)
    struct = p.get_structure(pdb_id, pdb_fn)
    sr = ShrakeRupley()
    sr.compute(struct, level="S")
    return struct.sasa

def get_metric_plot_variables(pdb_id, gt_pdb_fn, pred_pdb_fn, ranking_fn, chain_ids, intrfc_dist):
    ''' Get variables that we will calculate dockq/rmsd value against
        @Return idr_len, num_interface_resid, plddt, sasa, helix chain
    '''
    reverted = chain_ids[-1]
    chain_ids = chain_ids[:-1]

    load_and_select(intrfc_dist, gt_pdb_fn, pred_pdb_fn, chain_ids)
    cmd.super('native_R','pred_R')

    # len of ligand
    resid = read_residue_from_pdb(pred_pdb_fn)
    if reverted: # [ligand,receptor]
        len_ligand = len(resid[0])
    else:
        len_ligand = len(resid[1])

    # num of interface resid for ligand
    intrfc_resids = extract_residue_from_selection('pred_interface_L')
    len_intrfc_resid = len(intrfc_resids)

    # plddt
    fp = open(ranking_fn)
    rank = json.load(fp)
    plddt = rank["plddts"]["model_1_pred_0"]

    # comprehensive sasa: (sasa_r + sasa_l - sasa_super)/2
    #sasa = get_sasa(pdb_id, gt_pdb_fn) + get_sasa(pdb_id, pred_pdb_fn) - get_sasa(pdb_id, complex_fn)
    #sasa /= 2
    sasa = 0
    return (len_ligand, len_intrfc_resid, plddt, sasa)
