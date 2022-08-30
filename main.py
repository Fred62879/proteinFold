
import os
import utils
import pickle
import argparse
import numpy as np
import configargparse

from pymol import cmd
from parser import parse_args
from os.path import join, exists

def process_input(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''

    chain_start_resid_ids = utils.generate_fasta_from_pdb\
        (args.input_pdb_dir, args.linker_fasta_dir, args.pdb_ids, args.n_g)

    gt_chain_bd_resid_ids = utils.read_bd_resid_id_all \
        (args.input_pdb_dir, args.pdb_ids)

    with open(args.chain_start_resid_ids_fn, 'wb') as fp:
        pickle.dump(chain_start_resid_ids, fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'wb') as fp:
        pickle.dump(gt_chain_bd_resid_ids, fp)

    # get chain names and store locally
    chain_names = utils.get_chain_identifiers_all \
        (args.pdb_ids, args.input_pdb_dir)
    with open(args.chain_names_fn, 'wb') as fp:
        pickle.dump(chain_names, fp)


def process_output(args):

    with open(args.chain_names_fn, 'rb') as fp:
        chain_names = pickle.load(fp)

    with open(args.chain_start_resid_ids_fn, 'rb') as fp:
        chain_start_ids = pickle.load(fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'rb') as fp:
        gt_chain_bd_ids = pickle.load(fp)

    rmsds = []
    for pdb_id in args.pdb_ids[0:1]:
        print(pdb_id)
        dir = join(args.output_dir, pdb_id + '.fasta')
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            continue

        pred_fn = join(dir, args.pred_fn_str)
        gt_pdb_fn = join(args.input_pdb_dir, pdb_id+ '.pdb')
        pred_removed_linker_fn = join(dir, args.removed_linker_fn_str)

        utils.remove_linker(pred_fn, pred_removed_linker_fn, args.n_g,
                            chain_names[pdb_id], chain_start_ids[pdb_id],
                            gt_chain_bd_ids[pdb_id])

        # prune unknown aa from sequences
        if args.prune_unknown:
            assert(False)

        cur_rmsds = utils.calculate_rmsd\
            (gt_pdb_fn, pred_removed_linker_fn, chain_names[pdb_id],
             backbone=args.backbone, remove_hydrogen=args.remove_hydrogen)
        rmsds.append(cur_rmsds)

    print(args.rmsd_fn)
    np.savetxt(args.rmsd_fn, np.array(rmsds),delimiter=',',header=args.rmsd_names,comments='')
    cmd.quit()


if __name__ == "__main__":

    parser = configargparse.ArgumentParser()
    config = parse_args(parser)
    args = argparse.Namespace(**config)

    for operation in args.operations:
        if operation == 'process_input':
            process_input(args)
        elif operation == 'process_output':
            process_output(args)
