
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

def write_to_csv(data, out_fn):
    with open(out_fn, 'w', newline='') as fp:
        writer = csv.writer(fp)
        for row in data:
            writer.writerow(row)

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
    seq = []
    pdb_parser = PDBParser(QUIET = True)
    structure = pdb_parser.get_structure("model", fn)
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
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
    ''' find residue id of two boundary residues in each chain '''
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

def parse_pdb_ids(dir, suffix):
    # collect all unique pdb ids in a given directory
    pdb_ids = list(set(listdir(dir))) #next(os.walk(dir))[1]
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

#########
# metric
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
    best_model_id = rank['order'][0]
    plddt = rank['plddts'][best_model_id]

    # comprehensive sasa: (sasa_r + sasa_l - sasa_super)/2
    #sasa = get_sasa(pdb_id, gt_pdb_fn) + get_sasa(pdb_id, pred_pdb_fn) - get_sasa(pdb_id, complex_fn)
    #sasa /= 2
    sasa = 0
    return (len_ligand, len_intrfc_resid, plddt, sasa)

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
