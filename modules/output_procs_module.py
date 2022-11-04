
import os
import pickle
import numpy as np

from pymol import cmd
from os.path import join, exists
from biopandas.pdb import PandasPdb


class pipeline:
    def __init__(self, option, args):
        self.verbose = args.verbose

        self.n_g = args.n_g
        self.pdb_ids = args.pdb_ids
        self.input_pdb_dir = args.input_pdb_dir
        self.chain_ids_fn = args.chain_ids_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn

        self.pdb_download_command = args.pdb_download_command

        if option == 'init_output_procs':
            self._init_output_processing(args)

    def _init_output_processing(self, args):
        ''' pred_fn is filename of alphafold prediction
            if predict using fasta, the processing output is stored in aligned_fn
            otherwise stored in removed_linker_fn
        '''
        self.dockq = args.dockq
        self.backbone = args.backbone
        self.from_fasta = args.from_fasta
        self.remove_hydrogen = args.remove_hydrogen
        self.prune_and_renumber = args.prune_and_renumber
        self.interface_dist = args.interface_dist
        self.metric_col_names = args.metric_col_names

        self.fnat_path = args.fnat_path
        self.output_dir = args.output_dir
        self.metric_fn = args.metric_fn
        self.pred_fn_str = args.pred_fn_str
        self.ranking_fn_str = args.ranking_fn_str
        self.complex_fn_str = args.complex_fn_str
        self.aligned_fn_str = args.aligned_fn_str
        self.removed_linker_fn_str = args.removed_linker_fn_str

    def select_pdb_ids(config):
        #config['pdb_ids'] = ['1S5R']
        #config['pdb_ids'] = ['1C8O','1CF4','1D5S']
        #config['pdb_ids'] = ['2GGV','3H8K','5L0T']
        #config['pdb_ids'] = ['1JMT', '2AFF', '2DOH', '2GP9', '2ZQP', '3HK3', '3M51', '3MN7', '4LX3','6SBA','7BY2','7OY3']

        fn = config['pdb_ids_fn']
        excld_fn = config['pdb_exclude_fn']
        gpu_done_fn = config['pdb_gpu_done_fn']

        if exists(gpu_done_fn):
            pdb_ids = np.load(gpu_done_fn)
            pdb_ids.sort()
        elif exists(fn):
            pdb_ids = np.load(fn)
            if exists(excld_fn):
                exclude_ids = np.load(excld_fn)
                pdb_ids = np.array(list( set(pdb_ids) - set(exclude_ids) ))
            pdb_ids = np.sort(pdb_ids)
        else:
            all_pdb_ids = np.array(os.listdir(self.input_pdb_dir))
            pdb_ids = np.random.choice(all_pdb_ids, config['n_seqs'], replace=False)
            split_char = '_' if config['from_fasta'] and config['separate_fasta'] else '.'
            pdb_ids = np.array([id.split(split_char)[0] for id in pdb_ids])
            pdb_ids = np.sort(pdb_ids)
            np.save(fn, pdb_ids)

        config['pdb_ids'] = pdb_ids
        print(f'selected proteins {pdb_ids}')

    # chain_id is ordered as [receptor,ligand] which may differ from its native order
    def chain_id_orig_order(self, chain_ids):
        ligand_receptor_reverted = chain_ids[-1]
        chain_ids = chain_ids[:-1]
        if ligand_receptor_reverted:
            chain_ids = chain_ids[::-1]
        return chain_ids

    def remove_linker(self, pdb_id, pred_fn, pred_removed_linker_fn):
        ''' Remove linker from predicted pdb and rename chain
              identifier as per gt pdb
              i) if aa seq comes from gt pdb then also renumber
                 residue id as per gt pdb
             ii) if aa seq comes from fasta, then reorder residue id
                 such that after removal, residue id is contiguous (1,2,...,n)

            Note: pred pdb residue id is contiguous 1-indexed
              and since all chains are polyg linked, resid id
              is contiguous between chains (polyg linker inclusive)

            Input: chain_ids   identifier for all chain of current gt pdb
                   chain_start_ids id of first residue in each chain (from pred pdb)
                   gt_chain_bd_ids id of two boundary residue of each chain (from gt pdb)
        '''
        #if exists(pred_removed_linker_fn): return

        # chain names is changed during input processing to [receptor,ligand]
        # while linker removal needs original order of chain names
        chain_ids = self.chain_id_orig_order(self.chain_ids[pdb_id])
        n = len(chain_ids)

        chain_start_ids = self.chain_start_ids[pdb_id]
        gt_chain_bd_ids = self.gt_chain_bd_ids[pdb_id]

        residues = utils.read_residue_from_pdb(pred_fn)[0]
        for i in range(1, len(chain_start_ids)-1):
            id = chain_start_ids[i] - 1 # 1-based to 0-based
            assert(residues[id - self.n_g : id] == 'G' * self.n_g)

        ppdb = PandasPdb()
        _ = ppdb.read_pdb(pred_fn)
        df = ppdb.df['ATOM']

        acc_linker_length = 0
        linker_atom_ids = []
        atom_los, atom_his, offsets = [], [], []

        for i in range(n): # for each chain
            resid_lo = chain_start_ids[i]
            # resid_hi is residue id of first 'g' in pred
            resid_hi = chain_start_ids[i+1] - self.n_g

            # make sure gt and pred chain has same length
            gt_resid_lo, gt_resid_hi = gt_chain_bd_ids[i]
            if not self.from_fasta:
                # if pred based on fasta, sequence length might differ
                # needs additional processing after linker removal
                assert(gt_resid_hi - gt_resid_lo + 1 == resid_hi - resid_lo)

            # locate all atoms in pred pdb for current chain
            atom_ids = df[ (df['residue_number'] >= resid_lo) &
                           (df['residue_number'] < resid_hi)].index
                           #(df['chain_id'] == chain_ids[i]) ].index
            ''' Don't perform chain_id comparison here since linked pdb has
                one single chain. However, this is necessary for gt pdb where
                different chain may have residues with same id.
            '''

            atom_lo, atom_hi = atom_ids[0], atom_ids[-1]

            ''' reorder residue id and rename chain
                df.loc upper bound is inclusive
                atom_hi is id of last atom of current chain
            '''
            if self.from_fasta:
                offset = acc_linker_length
                acc_linker_length += self.n_g
            else:
                offset = resid_lo - gt_resid_lo

            atom_los.append(atom_lo);atom_his.append(atom_hi)
            offsets.append(offset)
            df.loc[atom_lo:atom_hi ,'chain_id'] = chain_ids[i]

            # mark polyg linker for removal
            if i != n - 2:
                linker_atom_lo = atom_hi+1
                linker_resid_hi = chain_start_ids[i+1]
                linker_atom_hi = df[ (df['residue_number'] < linker_resid_hi) ].index[-1]
                linker_atom_ids += list(np.arange(linker_atom_lo, linker_atom_hi+1))

        # re-number chain residue id
        for atom_lo, atom_hi, offset in zip(atom_los, atom_his, offsets):
            df.loc[atom_lo:atom_hi ,'residue_number'] -= offset

        # drop linker residues
        df.drop(linker_atom_ids, axis=0, inplace=True)
        ppdb.df['ATOM'] = df
        ppdb.to_pdb(pred_removed_linker_fn)

    def remove_extra_residue_and_renumber(self, pdb_id, aligned_fn, removed_linker_fn, gt_pdb_fn, pred_removed_linker_fn):
        chain_ids = self.chain_id_orig_order(self.chain_ids[pdb_id])
        seqs1 = utils.read_residue_from_pdb(removed_linker_fn)
        seqs2 = utils.read_residue_from_pdb(gt_pdb_fn)
        ranges = utils.find_prune_ranges_all_chains(seqs1, seqs2, chain_ids)
        utils.prune_renumber_seq_given_ranges(pred_removed_linker_fn, aligned_fn,
                                              chain_ids, ranges[0], self.gt_chain_bd_ids[pdb_id])

    def calculate_rmsd(self, pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn, verbose=False):
        ''' Calculate rmsd between gt and pred
              superimpose pred onto gt (with receptor only)
              chain_ids here is ordered as [receptor,ligand]
        '''
        rmsds = []
        chain_ids = self.chain_ids[pdb_id][:-1]
        utils.load_and_select \
            (self.interface_dist, gt_pdb_fn, pred_pdb_fn,
             chain_ids, backbone=self.backbone,
             remove_hydrogen=self.remove_hydrogen)

        # superimpose receptor chains and calculate rmsd for ligand
        utils.superimpose_receptors(complex_fn)
        rmsd = cmd.rms_cur('native_L','pred_L')
        rmsds.append(rmsd)

        # save two objects after superimposing receptor chain
        cmd.color('gray','native')
        cmd.color('red','pred')
        for obj in ['native','pred']:
            cmd.color('yellow', f'{obj}_interface_R')
            cmd.color('blue',f'{obj}_interface_L')

        # calculate rmsd for interface idr only
        rmsd = cmd.rms_cur('native_interface_L','pred_interface_L')
        rmsds.append(rmsd)

        rmsds = [round(rmsd,3) for rmsd in rmsds]
        if verbose: print(rmsds)
        return rmsds

    def process_output_for_one_pdb(self, pdb_id, chain_ids):
        print(f'\r\nProcessing output for {pdb_id}')
        dir = join(self.output_dir, pdb_id + '.fasta')
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            return

        pred_fn = join(dir, self.pred_fn_str)
        ranking_fn = join(dir, self.ranking_fn_str)
        gt_pdb_fn = join(self.input_pdb_dir, pdb_id + '.pdb')
        pred_removed_linker_fn = join(dir, self.removed_linker_fn_str)
        complex_fn = join(self.output_dir, pdb_id + '.fasta', self.complex_fn_str)

        self.remove_linker(pdb_id, pred_fn, pred_removed_linker_fn)

        # remove extra residues from fasta aa so that pred and gt have same aa sequence
        # also renumber pred residue id so that it's the same as per gt
        if self.from_fasta and self.prune_and_renumber:
            pred_aligned_fn = join(dir, self.aligned_fn_str)
            self.remove_extra_residue_and_renumber \
                (pdb_id, pred_aligned_fn, pred_removed_linker_fn, gt_pdb_fn,
                 pred_removed_linker_fn)

        if not self.from_fasta or not self.prune_and_renumber:
            pred_pdb_fn = pred_removed_linker_fn
        else: pred_pdb_fn = pred_aligned_fn

        if self.dockq:
            dockq_metrics = calc_metrics(pred_pdb_fn, gt_pdb_fn, self.fnat_path)
            (irms, Lrms, dockQ) = dockq_metrics['irms'], dockq_metrics['Lrms'], dockq_metrics['dockQ']
            #print(irms, Lrms, dockQ)
            cur_metrics = [pdb_id, irms, Lrms, dockQ]
        else:
            cur_metrics = self.calculate_rmsd \
                (pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn, self.verbose)
            cur_metrics.insert(0, pdb_id)

        variables = utils.get_metric_plot_variables(pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn,
                                                    ranking_fn, chain_ids, self.interface_dist)
        cur_metrics = np.append(cur_metrics, variables)
        return cur_metrics

    def process_output(self):
        with open(self.chain_ids_fn, 'rb') as fp:
            self.chain_ids = pickle.load(fp)

        with open(self.chain_start_resid_ids_fn, 'rb') as fp:
            self.chain_start_ids = pickle.load(fp)

        with open(self.gt_chain_bd_resid_ids_fn, 'rb') as fp:
            self.gt_chain_bd_ids = pickle.load(fp)

        failed_pdbs = []
        metrics = [self.metric_col_names]
        for pdb_id in self.pdb_ids:
            chain_ids = self.chain_ids[pdb_id]
            #cur_metrics = self.process_output_for_one_pdb(pdb_id, chain_ids)
            #metrics.append(cur_metrics)

            try:
                cur_metrics = self.process_output_for_one_pdb(pdb_id, chain_ids)
            except Exception as e:
                failed_pdbs.append(pdb_id)
                print('! error', pdb_id, e)
            else:
                metrics.append(cur_metrics)

        utils.write_to_csv(metrics, self.metric_fn)

        pdb_ids, metric_names, data = utils.parse_csv(self.metric_fn)
        utils.plot_scatter(self.output_dir, data, pdb_ids)
        print('pdbs failed: ', failed_pdbs)
