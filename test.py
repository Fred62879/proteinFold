
import os
import utils

from pymol import cmd
from os.path import join
from biopandas.pdb import PandasPdb


def run(input_dir, output_dir, gt_pdb_fn, pred_pdb_fn, chains=None):
    pdb_id = gt_pdb_fn.split('.')[0]
    print(input_dir, pdb_id)

    cmd.load(join(input_dir, gt_pdb_fn), 'native')
    #cmd.png(join(output_dir, 'native.png'))

    chains = split_chains()
    print(chains)

    #cmd.remove('to_remove')

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, align_chains):
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    cmd.select('pred_A', 'pred and chain A')
    cmd.select('pred_B', 'pred and chain B')
    cmd.select('native_A', 'native and chain A')
    cmd.select('native_B', 'native and chain B')
    B_super = cmd.super('pred_B','native_B')
    B_align = cmd.align('pred_B','native_B')
    print(B_super


def main(input_dir, output_dir):
    fns = os.listdir(input_dir)
    print(fns)
    for gt_pdb_fn in fns[0:5]:
        if 'pdb' not in gt_pdb_fn: continue

        pred_pdb_fn = ''
        run(input_dir, output_dir, gt_pdb_fn, pred_pdb_fn)

    cmd.quit()

if __name__ == "__main__":

    input_dir = '/media/fred/Local Disk/Projects/bioinfo/data/input/idr_84/pdbs'
    #input_dir = '/media/fred/Local Disk/Projects/bioinfo/data/input/peptide/pdbs'
    output_dir = '/media/fred/Local Disk/Projects/bioinfo/data/output'

    #main(input_dir, output_dir)
    #print(utils.split_chains_wo_env(join(input_dir, '1YCQ.pdb')))
