#import sys
import pickle

from pymol import cmd
#sys.path.insert(0,'/media/fred/Local Disk/bioinfo/code/pipeline')
#from utils import load_and_select

def load_select_super(pdb_id, interface_dist):
    # data_dir = '/media/fred/working_drive/projects/bioinfo/data/'
    data_dir = '/scratch/projects/bioinfo/data'
    pdb_fn = f'{data_dir}/input/ds1/pdbs/{pdb_id}.pdb'
    pred_fn = f'{data_dir}/output/ds1_af_full/poly_g_20_fasta/{pdb_id}.fasta/ranked_0_removed_linker_aligned.pdb'
    cmd.load(pdb_fn, 'native')
    cmd.load(pred_fn, 'pred')

    with open(f'{data_dir}/input/ds1/poly_g_20_fasta/chain_names0.pkl','rb') as fp:
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

def load_and_select(interface_dist, gt_pdb_fn, pred_pdb_fn, ordered_chain_ids, backbone=False, remove_hydrogen=False):
    # Load gt and pred pdb files and select receptor, ligand, and interface
    cmd.delete('all')
    cmd.load(gt_pdb_fn, 'native')
    cmd.remove('native and solvent')
    cmd.remove('native and hydrogens')
    cmd.load(pred_pdb_fn, 'pred')

    # remove hydrogens (presented in af prediction)
    if remove_hydrogen:
        cmd.remove('hydrogens')

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


ds = 'ds2'
pdb_id = '1A2X'
interface_dist = 5

# data_dir = '/media/fred/working_drive/projects/bioinfo/data/'
# pdb_fn = f'/media/fred/working_drive/projects/bioinfo/data/input/ds2/pdbs/{pdb_id}.pdb'
# pred_fn = f'/media/fred/working_drive/projects/bioinfo/data/output/ds2_af_full/poly_g_20/{pdb_id}.fasta/ranked_0_aligned.pdb'

data_dir = '/scratch/projects/bioinfo/data'
pdb_fn = f'/scratch/projects/bioinfo/data/input/colabfold_evaluate/pdbs/{pdb_id}.pdb'
# pred_fn = f'/scratch/projects/bioinfo/data/output/multimer_ptm/{pdb_id}_unrelaxed_rank_1_model_4.pdb'
pred_fn = f'/scratch/projects/bioinfo/data/output/ptm_poly_g/{pdb_id}__unknown_description__unrelaxed_rank_1_model_4.pdb'

with open(f'{data_dir}/input/{ds}/ordered_chain_ids.pkl','rb') as fp:
    ids = pickle.load(fp)

ids = ids[pdb_id][:2] # [receptor,ligand]
load_and_select(interface_dist, pdb_fn, pred_fn, ids, backbone=True)

# run /media/fred/working_drive/projects/bioinfo/code/pipeline/utils/pymol_demo.py
# run /scratch/projects/bioinfo/code/pipeline/utils/pymol_demo.py
