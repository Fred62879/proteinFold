
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
        self.failed_pdbs_fn = args.failed_pdbs_fn
        self.pdb_exclude_fn = args.pdb_exclude_fn
        self.pdb_gpu_done_fn = args.pdb_gpu_done_fn
        self.orig_chain_ids_fn = args.orig_chain_ids_fn
        self.ordered_chain_ids_fn = args.ordered_chain_ids_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn
        self.second_structs_resid_ids_fn = args.second_structs_resid_ids_fn

        self.pruned_fn_str = args.pruned_fn_str
        self.ranking_fn_str = args.ranking_fn_str
        self.aligned_fn_str = args.aligned_fn_str
        self.pred_pdb_fn_str = args.pred_pdb_fn_str
        self.superimposed_fn_str = args.superimposed_fn_str
        self.removed_linker_fn_str = args.removed_linker_fn_str

        with open(self.orig_chain_ids_fn, 'rb') as fp:
            self.orig_chain_ids = pickle.load(fp)
        with open(self.ordered_chain_ids_fn, 'rb') as fp:
            self.ordered_chain_ids = pickle.load(fp)

        if exists(self.second_structs_resid_ids_fn):
            with open(self.second_structs_resid_ids_fn, 'rb') as fp:
                self.second_structs_resid_ids = pickle.load(fp)
        else:
            self.second_structs_resid_ids = None

        if self.strategy == 'poly_g_link':
            with open(self.chain_start_resid_ids_fn, 'rb') as fp:
                self.chain_start_ids = pickle.load(fp)
        if self.strategy == 'poly_g_link' or self.prune_extra_residues or self.renumber_residues:
            with open(self.gt_chain_bd_resid_ids_fn, 'rb') as fp:
                self.gt_chain_bd_ids = pickle.load(fp)

    def process_output(self):
        pdb_ids = self.select_pdb_ids()
        failed_pdbs, metrics = {}, [self.metric_col_names]

        if exists(self.metric_fn):
            _, metric_names, metrics_done = utils.parse_csv(self.metric_fn)
        else: metrics_done = None

        for pdb_id in pdb_ids:
            cur_pdb_dir = join(self.output_dir, pdb_id + '.fasta')
            gt_pdb_fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            pred_pdb_fn = join(cur_pdb_dir, self.pred_pdb_fn_str + ".pdb")

            pred_pdb_fn = self.process_predicted_pdb(pdb_id, gt_pdb_fn, pred_pdb_fn)
            cache = None if metrics_done is None or pdb_id not in metrics_done else metrics_done[pdb_id]
            select_residue_ids = None if self.second_structs_resid_ids is None \
                else self.second_structs_resid_ids[pdb_id]

            cur_metrics = self.calculate_metrics(
                cur_pdb_dir, pdb_id, gt_pdb_fn, pred_pdb_fn, cache,
                selected_residue_ids=select_residue_ids)
            metrics.append(cur_metrics)

            '''
            try:
                pred_pdb_fn = self.process_predicted_pdb(pdb_id, gt_pdb_fn, pred_pdb_fn)
                if metrics_done is not None:
                    cache = None if pdb_id not in metrics_done else metrics_done[pdb_id]
                    cur_metrics = self.calculate_metrics(cur_pdb_dir, pdb_id, gt_pdb_fn, pred_pdb_fn, cache)
                else:
                    select_residue_ids = None if self.second_struct_ids is None \
                        else self.second_struct_ids[pdb_id]
                    cur_metrics = self.calculate_metrics(cur_pdb_dir, pdb_id, gt_pdb_fn, pred_pdb_fn,
                                                         selected_residue_ids=select_residue_ids)
            except Exception as e:
                failed_pdbs[pdb_id] = f'{e}'
                print(f'ERROR: {pdb_id} {e}')
            else:
                metrics.append(cur_metrics)
            '''

        utils.write_to_csv(metrics, self.metric_fn)
        pdb_ids, metric_names, data = utils.parse_csv(self.metric_fn)
        utils.plot_scatter(self.output_dir, data, pdb_ids)
        #print('pdbs failed: ', failed_pdbs)
        with open(self.failed_pdbs_fn, 'wb') as fp:
            pickle.dump(failed_pdbs, fp)

    ##########
    # helpers
    def select_pdb_ids(self):
        excld_fn = self.pdb_exclude_fn
        gpu_done_fn = self.pdb_gpu_done_fn

        pdb_ids = utils.parse_pdb_ids(self.output_dir, '.fasta')
        '''
        if exists(gpu_done_fn):
            pdb_ids = np.load(gpu_done_fn)
        else:
            pdb_ids = utils.parse_pdb_ids(self.output_dir, '.fasta')
            #pdb_ids = np.random.choice(all_pdb_ids, self.n_seqs, replace=False)
        '''
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
            chain_ids = self.orig_chain_ids[pdb_id]
            removed_linker_fn = join(dir, f"{self.removed_linker_fn_str}_{chain_ids[0]}_{chain_ids[1]}.pdb")
            outils.remove_linker(pdb_id, pred_pdb_fn, removed_linker_fn,
                                 chain_ids, self.chain_start_ids[pdb_id],
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

    def calculate_metrics(self, dir, pdb_id, gt_pdb_fn, pred_pdb_fn, metrics_done=None, selected_residue_ids=None):
        """ Calculate metrics
            @Param:
               residue_ids: if not None, designate residues used for metric calculation
        """

        # whether the order of chains in gt pdb file is receptor-ligand. if not, set reverted as True
        reverted = self.ordered_chain_ids[pdb_id] != self.orig_chain_ids[pdb_id]

        if self.dockq:
            if metrics_done is not None and 'irms' in metrics_done and 'Lrms' in metrics_done and 'dockQ' in metrics_done:
                # (irms, Lrms, dockQ) = metrics_done['irms'], metrics_done['Lrms'], metrics_done['dockQ']
                cur_metrics = metrics_done
            else:
                cur_metrics = [pdb_id]

                if selected_residue_ids is not None:
                    for second_struct_ids in selected_residue_ids: # each 2nd structure
                        selected_ligand_resid_ids = second_struct_ids[1] if reverted else second_struct_ids[0]

                        dockq_metrics = calc_metrics(pred_pdb_fn, gt_pdb_fn, self.fnat_path,
                                                     selected_residue_ids=set(selected_ligand_resid_ids))

                        # metrics based off ligand residues within current 2nd structure
                        cur_metrics += [dockq_metrics['irms'], dockq_metrics['Lrms'], dockq_metrics['dockQ']]

                # metircs without considering 2nd structures
                dockq_metrics = calc_metrics(pred_pdb_fn, gt_pdb_fn, self.fnat_path)
                (irms, Lrms, dockQ) = dockq_metrics['irms'], dockq_metrics['Lrms'], dockq_metrics['dockQ']
                cur_metrics += [irms, Lrms, dockQ]
        else:
            if metrics_done is not None and 'rmsd' in metrics_done:
                cur_metrics = [metrics_done['rmsd']]
            else:
                cur_metrics = outils.calculate_rmsd \
                (pdb_id, gt_pdb_fn, pred_pdb_fn, self.interface_dist,
                 self.order_chain_ids[pdb_id], self.remove_backbone, self.remove_hydrogen)
            cur_metrics.insert(0, pdb_id)

        ranking_fn = join(dir, self.ranking_fn_str)
        superimposed_fn = join(dir, self.superimposed_fn_str)

        variables = utils.get_metric_plot_variables \
            (pdb_id, gt_pdb_fn, pred_pdb_fn, ranking_fn, superimposed_fn,
             self.ordered_chain_ids[pdb_id], self.interface_dist,
             reverted=reverted, metrics_done=metrics_done)

        cur_metrics = np.append(cur_metrics, variables)
        return cur_metrics
