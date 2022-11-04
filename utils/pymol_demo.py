#import sys
import pickle

from pymol import cmd
#sys.path.insert(0,'/media/fred/Local Disk/bioinfo/code/pipeline')
#from utils import load_and_select

def load_select_super(pdb_id, interface_dist):
    pdb_fn = f'/media/fred/Local Disk/projects/bioinfo/data/input/ds1/pdbs/{pdb_id}.pdb'
    pred_fn = f'/media/fred/Local Disk/projects/bioinfo/data/output/ds1_af_full/poly_g_20_fasta/{pdb_id}.fasta/ranked_0_removed_linker_aligned.pdb'
    cmd.load(pdb_fn, 'native')
    cmd.load(pred_fn, 'pred')
    with open('/media/fred/Local Disk/projects/bioinfo/data/input/ds1/poly_g_20_fasta/chain_names0.pkl','rb') as fp:
        ids = pickle.load(fp)
    id = ids[pdb_id] # [receptor,ligand]

    cmd.select('nativeA', f'native and chain {id[0]}')
    cmd.select('nativeB', f'native and chain {id[1]}')
    cmd.select('predA', f'pred and chain {id[0]}')
    cmd.select('predB', f'pred and chain {id[1]}')
    cmd.color('orange','nativeA')
    cmd.color('blue','nativeB')
    cmd.color('yellow','predA')
    cmd.color('gray','predB')

    cmd.super('nativeA','predA')

#load_select_super('4MMT')

def set_chain_selector(receptor_chain_id, ligand_chain_id, backbone=False):
    ligand_chain_selector = f'chain {ligand_chain_id}'
    receptor_chain_selector = f'chain {receptor_chain_id}'

    if backbone:
        ligand_chain_selector += ' and backbone'
        receptor_chain_selector += ' and backbone'

    return receptor_chain_selector, ligand_chain_selector

def load_and_select(interface_dist, gt_pdb_fn, pred_pdb_fn, chain_ids, backbone=False, remove_hydrogen=False):
    # Load gt and pred pdb files and select receptor, ligand, and interface
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.remove('native and solvent')
    cmd.remove('native and hydrogens')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

    # only deal with dimers
    assert(len(chain_ids) == 2)

    # first chain is receptor
    receptor_chain_selector, ligand_chain_selector = set_chain_selector \
        (chain_ids[0], chain_ids[1], backbone)

    for obj in ['native','pred']:
        # select ligand and receptor based on initial assumption
        cmd.select(f'{obj}_R', f'{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_L', f'{obj} and {ligand_chain_selector}')

        # select receptor interface atoms
        cmd.select(f'{obj}_R_interface_init', f'byres {obj}_R within {interface_dist} of {obj}_L')
        cmd.select(f'interface_{obj}', f'{obj}_R_interface_init + byres {obj}_L within {interface_dist} of {obj}_R_interface_init')

        cmd.select(f'{obj}_interface_R', f'interface_{obj} and {receptor_chain_selector}')
        cmd.select(f'{obj}_interface_L', f'interface_{obj} and {ligand_chain_selector}')

    cmd.color('gray','native_R')
    cmd.color('blue','native_L')
    cmd.color('brown','pred_R')
    cmd.color('purple','pred_L')

    cmd.color('orange', 'native_interface_R')
    cmd.color('green','native_interface_L')
    cmd.color('orange', 'pred_interface_R')
    cmd.color('green','pred_interface_L')

    cmd.super('native_R','pred_R')
    LRMS = cmd.rms_cur('native_L','pred_L')
    iRMS = cmd.align('native_interface_R','pred_interface_R')
    print(LRMS, iRMS)

pdb_id = '5L0T'
interface_dist = 10
pdb_fn = f'/media/fred/Local Disk/projects/bioinfo/data/input/ds1/pdbs/{pdb_id}.pdb'
pred_fn = f'/media/fred/Local Disk/projects/bioinfo/data/output/ds1_af_full/poly_g_20_fasta/{pdb_id}.fasta/ranked_0_removed_linker_aligned.pdb'
with open('/media/fred/Local Disk/projects/bioinfo/data/input/ds1/poly_g_20_fasta/chain_ids0.pkl','rb') as fp:
    ids = pickle.load(fp)
ids = ids[pdb_id][:2] # [receptor,ligand]
load_and_select(interface_dist, pdb_fn, pred_fn, ids, backbone=True)

#run /media/fred/Local Disk/projects/bioinfo/code/pipeline/utils/pymol_demo.py
