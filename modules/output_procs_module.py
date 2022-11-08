
import os
import pickle
import numpy as np
import utils.common as utils
import utils.output_procs as outils

from pymol import cmd
from os.path import join, exists
from utils.dockq import calc_metrics


class pipeline:
    def __init__(self, option, args):

        ''' pred_fn is filename of alphafold prediction
            if predict using fasta, the processing output is stored in aligned_fn
            otherwise stored in removed_linker_fn
        '''
        self.strategy = args.strategy

        self.dockq = args.dockq
        self.verbose = args.verbose
        self.backbone = args.backbone
        self.remove_hydrogen = args.remove_hydrogen
        self.renumber_residues = args.renumber_residues
        self.prune_extra_residues = args.prune_extra_residues
        self.generate_fasta_from_pdb = args.generate_fasta_from_pdb

        self.n_g = args.n_g
        self.interface_dist = args.interface_dist
        self.metric_col_names = args.metric_col_names

        self.fnat_path = args.fnat_path
        self.metric_fn = args.metric_fn
        self.output_dir = args.output_dir
        self.input_pdb_dir = args.input_pdb_dir
        self.pdb_exclude_fn = args.pdb_exclude_fn
        self.pdb_gpu_done_fn = args.pdb_gpu_done_fn
        self.orig_chain_ids_fn = args.orig_chain_ids_fn
        self.ordered_chain_ids_fn = args.ordered_chain_ids_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn

        self.ranking_fn_str = args.ranking_fn_str
        self.pruned_fn_str = args.pruned_fn_str
        self.aligned_fn_str = args.aligned_fn_str
        self.pred_pdb_fn_str = args.pred_pdb_fn_str
        self.removed_linker_fn_str = args.removed_linker_fn_str

        with open(self.orig_chain_ids_fn, 'rb') as fp:
            self.orig_chain_ids = pickle.load(fp)
        with open(self.ordered_chain_ids_fn, 'rb') as fp:
            self.ordered_chain_ids = pickle.load(fp)

        if self.strategy == 'poly_g_link':
            with open(self.chain_start_resid_ids_fn, 'rb') as fp:
                self.chain_start_ids = pickle.load(fp)
        if self.strategy == 'poly_g_link' or self.prune_extra_residues or self.renumber_residues:
            with open(self.gt_chain_bd_resid_ids_fn, 'rb') as fp:
                self.gt_chain_bd_ids = pickle.load(fp)

    def process_output(self):
        pdb_ids = self.select_pdb_ids()
        failed_pdbs, metrics = [], [self.metric_col_names]

        for pdb_id in pdb_ids:
            cur_pdb_dir = join(self.output_dir, pdb_id + '.fasta')
            gt_pdb_fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            pred_pdb_fn = join(cur_pdb_dir, self.pred_pdb_fn_str)

            #pred_pdb_fn = self.process_predicted_pdb(pdb_id, gt_pdb_fn, pred_pdb_fn)
            #cur_metrics = self.calculate_metrics(cur_pdb_dir, pdb_id, gt_pdb_fn, pred_pdb_fn)
            try:
                pred_pdb_fn = self.process_predicted_pdb(pdb_id, gt_pdb_fn, pred_pdb_fn)
                cur_metrics = self.calculate_metrics(cur_pdb_dir, pdb_id, gt_pdb_fn, pred_pdb_fn)
            except Exception as e:
                failed_pdbs.append(pdb_id)
                print(f'ERROR: {pdb_id} {e}')
            else: metrics.append(cur_metrics)

        utils.write_to_csv(metrics, self.metric_fn)
        pdb_ids, metric_names, data = utils.parse_csv(self.metric_fn)
        utils.plot_scatter(self.output_dir, data, pdb_ids)
        print('pdbs failed: ', failed_pdbs)

    ##########
    # helpers
    def select_pdb_ids(self):
        #pdb_ids = ['1S5R']
        excld_fn = self.pdb_exclude_fn
        gpu_done_fn = self.pdb_gpu_done_fn

        if exists(gpu_done_fn):
            pdb_ids = np.load(gpu_done_fn)
        else:
            pdb_ids = utils.parse_pdb_ids(self.output_dir, '.fasta')
            #pdb_ids = np.random.choice(all_pdb_ids, self.n_seqs, replace=False)

        if exists(excld_fn):
            exclude_ids = np.load(excld_fn)
            pdb_ids = np.array(list( set(pdb_ids) - set(exclude_ids) ))

        pdb_ids = np.sort(pdb_ids)
        print(f'selected {len(pdb_ids)} proteins')
        return pdb_ids

    def process_predicted_pdb(self, pdb_id, gt_pdb_fn, pred_pdb_fn):
        print(f'processing prediction for {pdb_id}')
        dir = join(self.output_dir, pdb_id + '.fasta')
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            return

        if self.strategy == 'poly_g_link':
            #print('  removing linker')
            removed_linker_fn = join(dir, self.removed_linker_fn_str)
            outils.remove_linker(pdb_id, pred_pdb_fn, removed_linker_fn,
                                 self.orig_chain_ids[pdb_id], self.chain_start_ids[pdb_id],
                                 self.gt_chain_bd_ids[pdb_id], self.n_g, self.generate_fasta_from_pdb)
            pred_pdb_fn = removed_linker_fn

        utils.assert_equal_num_chains(gt_pdb_fn, pred_pdb_fn)

        if self.prune_extra_residues:
            #print('  pruning extra residues')
            pruned_fn = join(dir, self.pruned_fn_str)
            outils.remove_extra_residue_predicted_pdb \
                (pred_pdb_fn, pruned_fn, gt_pdb_fn,
                 self.orig_chain_ids[pdb_id], self.gt_chain_bd_ids[pdb_id])
            pred_pdb_fn = pruned_fn

        if self.renumber_residues:
            #print('  renumbering residues')
            aligned_fn = join(dir, self.aligned_fn_str)
            outils.renumber_residue_perdicted_pdb \
                (pred_pdb_fn, aligned_fn, self.orig_chain_ids[pdb_id], self.gt_chain_bd_ids[pdb_id])
            pred_pdb_fn = aligned_fn

        utils.assert_seq_equality(gt_pdb_fn, pred_pdb_fn)
        return pred_pdb_fn

    def calculate_metrics(self, dir, pdb_id, gt_pdb_fn, pred_pdb_fn):
        if self.dockq:
            dockq_metrics = calc_metrics(pred_pdb_fn, gt_pdb_fn, self.fnat_path)
            (irms, Lrms, dockQ) = dockq_metrics['irms'], dockq_metrics['Lrms'], dockq_metrics['dockQ']
            cur_metrics = [pdb_id, irms, Lrms, dockQ]
        else:
            cur_metrics = outils.calculate_rmsd \
                (pdb_id, gt_pdb_fn, pred_pdb_fn, self.interface_dist,
                 self.order_chain_ids[pdb_id], self.remove_backbone, self.remove_hydrogen)
            cur_metrics.insert(0, pdb_id)

        ranking_fn = join(dir, self.ranking_fn_str)
        variables = utils.get_metric_plot_variables \
            (pdb_id, gt_pdb_fn, pred_pdb_fn, ranking_fn,
             self.ordered_chain_ids[pdb_id], self.interface_dist)

        cur_metrics = np.append(cur_metrics, variables)
        return cur_metrics
