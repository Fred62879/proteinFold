
import os
import utils
import pickle
import argparse
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
    #if args.remove_x:


    chain_start_resid_ids = utils.generate_fasta_from_pdb\
        (args.input_pdb_dir, args.linker_fasta_dir, args.pdb_ids, args.n_g,
         remove_x=args.remove_x)

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
    ''' remove linker from predicted pdb and reorder
          residue number as per gt pdb'
        calculate rmsd and export to csv file
    '''

    with open(args.chain_names_fn, 'rb') as fp:
        chain_names = pickle.load(fp)

    with open(args.chain_start_resid_ids_fn, 'rb') as fp:
        chain_start_ids = pickle.load(fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'rb') as fp:
        gt_chain_bd_ids = pickle.load(fp)

    rmsds = [args.rmsd_names]
    for pdb_id in args.pdb_ids:
        print()
        print(pdb_id, chain_names[pdb_id], gt_chain_bd_ids[pdb_id])

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

        complex_fn = join(args.output_dir, pdb_id + '.fasta', args.complex_fn_str)
        cur_rmsds = utils.calculate_rmsd\
            (gt_pdb_fn, pred_removed_linker_fn, complex_fn,
             chain_names[pdb_id], backbone=args.backbone,
             remove_hydrogen=args.remove_hydrogen)
        cur_rmsds.insert(0, pdb_id)
        rmsds.append(cur_rmsds)

    utils.write_to_csv(rmsds, args.rmsd_fn)
    cmd.quit()

def assert_fasta(args):
    ''' check if aa sequence from processed prediction pdb match
          sequence from gt pdb
    '''
    for id in args.pdb_ids:
        pdb_fn = join(args.input_pdb_dir, id + '.pdb')
        pred_fn = join(args.output_dir, id + '.fasta', args.pred_fn)
        seq_1 = utils.read_residue_from_pdb(pdb_fn)
        seq_2 = utils.read_residue_from_pdb(pred_fn)
        #print(seq_1)
        #print(seq_2)
        assert(seq_1 == seq_2)

if __name__ == "__main__":

    parser = configargparse.ArgumentParser()
    config = parse_args(parser)
    args = argparse.Namespace(**config)

    for operation in args.operations:
        if operation == 'process_input':
            process_input(args)
        elif operation == 'process_output':
            process_output(args)
        elif operation == 'assert_fasta':
            assert_fasta(args)
