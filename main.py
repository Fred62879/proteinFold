
import os
import utils
import pickle
import argparse
import configargparse

from pymol import cmd
from parser import parse_args
from os.path import join, exists


def process_output(args):

    with open(args.chain_names_fn, 'rb') as fp:
        chain_names = pickle.load(fp)

    with open(args.chain_start_resid_ids_fn, 'rb') as fp:
        chain_start_ids = pickle.load(fp)

    with open(args.gt_chain_bd_resid_ids_fn, 'rb') as fp:
        gt_chain_bd_ids = pickle.load(fp)

    for pdb_id in args.pdb_ids[0:1]:
        print()
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
            align_option = 'pruned_X'
            seqs_pred = utils.read_fasta_from_pdb(pred_removed_linker_fn)
            seqs_gt = utils.read_fasta_from_pdb(gt_pdb_fn)
            print(seqs_pred)
            print(seqs_gt)
            prune_ranges = utils.find_prune_ranges_all_chains(seqs_pred, seqs_gt, prune_X=True)

            #print(prune_ranges)
            gt_removed_X_fn = gt_pdb_fn[:-4] + '_' + align_option + '.pdb'
            pred_removed_X_fn = pred_fn[:-4] + '_' + align_option + '.pdb'
            utils.align_chains(gt_pdb_fn, gt_removed_X_fn, pred_removed_linker_fn, pred_removed_X_fn, prune_ranges)

        utils.calculate_rmsd(gt_pdb_fn, pred_removed_linker_fn, chain_names[pdb_id])

    #utils.peptide_metric(args.pdb_ids, args.input_pdb_dir, args.output_dir, args.pred_fn, chain_names)

    cmd.quit()

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

def assert_fasta(args):
    for id in args.pdb_ids:
        pdb_fn = join(args.input_pdb_dir, id + '.pdb')
        fasta_fn = join(args.source_fasta_dir, id + '.pdb')
        seq_1 = utils.read_residue_from_pdb(pdb_fn)
        seq_2 = utils.read_residue_from_pdb(fasta_fn)
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
