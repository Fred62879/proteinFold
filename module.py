
import os
import utils
import pickle
import numpy as np

from pymol import cmd
from functools import reduce
from dockq import calc_metrics
from os.path import join, exists
from biopandas.pdb import PandasPdb

class pipeline:
    def __init__(self, option, args):
        self.verbose = args.verbose

        self.n_g = args.n_g
        self.pdb_ids = args.pdb_ids
        self.input_pdb_dir = args.input_pdb_dir
        self.chain_names_fn = args.chain_names_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn

        if option == 'init_atom_prune':
            self._init_atom_prune(args)
        elif option == 'init_input_procs':
            self._init_input_processing(args)
        elif option == 'init_input_procs_fasta':
            self._init_input_processing_from_fasta(args)
        elif option == 'init_output_procs':
            self._init_output_processing(args)
        elif option == 'init_atom_locating':
            self._init_atom_locating(args)

    ###############
    # init methods
    def _init_input_processing(self, args):
        self.pdb_ids_fn = args.pdb_ids_fn
        self.linked_fasta_dir = args.linked_fasta_dir

    def _init_input_processing_from_fasta(self, args):
        self.source_fasta_dir = args.source_fasta_dir
        self.linked_fasta_dir = args.linked_fasta_dir

    def _init_output_processing(self, args):
        ''' pred_fn is filename of alphafold prediction
            if predict using fasta, the processing output is stored in aligned_fn
            otherwise stored in removed_linker_fn
        '''
        self.output_dir = args.output_dir
        self.rmsd_names = args.rmsd_names
        self.interface_dist = args.interface_dist
        self.fnat_path = args.code_dir

        self.dockq = args.dockq
        self.backbone = args.backbone
        self.from_fasta = args.from_fasta
        self.remove_hydrogen = args.remove_hydrogen
        self.prune_and_renumber = args.prune_and_renumber

        self.rmsd_fn = args.rmsd_fn
        self.pred_fn_str = args.pred_fn_str
        self.complex_fn_str = args.complex_fn_str
        self.aligned_fn_str = args.aligned_fn_str
        self.removed_linker_fn_str = args.removed_linker_fn_str

    def _init_atom_prune(self, args):
        self.prune_pdb_ids = args.prune_pdb_ids
        self.atom_prune_ranges = args.atom_prune_ranges

    def _init_atom_locating(self, args):
        self.pdb_ids = args.pdb_ids
        self.input_pdb_dir = args.input_pdb_dir

        self.from_fasta = args.from_fasta
        self.output_dir = args.output_dir
        self.pred_fn_str = args.pred_fn_str
        self.aligned_fn_str = args.aligned_fn_str

    #######################
    # processing functions
    def process_input(self):
        chain_start_resid_ids = self.generate_fasta_from_pdb()
        gt_chain_bd_resid_ids = self.read_bd_resid_id_all()

        with open(self.chain_start_resid_ids_fn, 'wb') as fp:
            pickle.dump(chain_start_resid_ids, fp)

        with open(self.gt_chain_bd_resid_ids_fn, 'wb') as fp:
            pickle.dump(gt_chain_bd_resid_ids, fp)

        # get chain names and store locally
        chain_names = self.read_chain_names_all()
        with open(self.chain_names_fn, 'wb') as fp:
            pickle.dump(chain_names, fp)

    def process_input_from_fasta(self):
        groups = self.parse_to_groups()
        chain_start_resid_ids = self.poly_g_link_all(groups)
        gt_chain_bd_resid_ids = self.read_bd_resid_id_all()

        with open(self.chain_start_resid_ids_fn, 'wb') as fp:
            pickle.dump(chain_start_resid_ids, fp)

        with open(self.gt_chain_bd_resid_ids_fn, 'wb') as fp:
            pickle.dump(gt_chain_bd_resid_ids, fp)

        # get chain names and store locally
        chain_names = self.read_chain_names_all()
        with open(self.chain_names_fn, 'wb') as fp:
            pickle.dump(chain_names, fp)

    def locate_extra_atoms(self):
        for pdb_id in self.pdb_ids:
            print(f'{pdb_id} extra atoms: ')
            pdb_fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            dir = join(self.output_dir, pdb_id + '.fasta')
            if self.from_fasta:
                pred_fn = join(dir, self.aligned_fn_str)
            else:
                pred_fn = join(dir, self.removed_linker_fn_str)
            ret = utils.find_residue_diff_in_atom_counts(pdb_fn, pred_fn)
            print(ret)

    def prune_extra_atoms(self):
        ''' prune_pdb_ids and prune_ranges are hardcoded in parser.py '''
        for pdb_id, rnge in zip(self.prune_pdb_ids, self.atom_prune_ranges):
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            utils.prune_pdb_atoms(fn, rnge)

    ''' Below are output processing methods '''
    def process_output_for_one_pdb(self, pdb_id):
        print(f'\r\nProcessing output for {pdb_id}')
        dir = join(self.output_dir, pdb_id + '.fasta')
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            return

        pred_fn = join(dir, self.pred_fn_str)
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
            out_fn = pred_removed_linker_fn
        else: out_fn = pred_aligned_fn

        if self.dockq:
            dockq_metrics = calc_metrics(out_fn, gt_pdb_fn, self.fnat_path)
            (irms, Lrms, dockQ) = dockq_metrics['irms'], dockq_metrics['Lrms'], dockq_metrics['dockQ']
            #print(irms, Lrms, dockQ)
            cur_metrics = [pdb_id, irms, Lrms, dockQ]
        else:
            cur_metrics = self.calculate_rmsd \
                (pdb_id, gt_pdb_fn, out_fn, complex_fn, self.verbose)
            cur_metrics.insert(0, pdb_id)
        return cur_metrics

    def process_output(self):
        with open(self.chain_names_fn, 'rb') as fp:
            self.chain_names = pickle.load(fp)

        with open(self.chain_start_resid_ids_fn, 'rb') as fp:
            self.chain_start_ids = pickle.load(fp)

        with open(self.gt_chain_bd_resid_ids_fn, 'rb') as fp:
            self.gt_chain_bd_ids = pickle.load(fp)

        failed_pdbs = []
        metrics = [self.rmsd_names]
        for pdb_id in self.pdb_ids:
            #cur_metrics = self.process_output_for_one_pdb(pdb_id)
            try:
                cur_metrics = self.process_output_for_one_pdb(pdb_id)
            except Exception as e:
                failed_pdbs.append(pdb_id)
                print(pdb_id, e)
            else:
                metrics.append(cur_metrics)

        utils.write_to_csv(metrics, self.rmsd_fn)
        print('pdbs failed: ', failed_pdbs)

        if not self.dockq:
            cmd.quit()

    def parse_to_groups(self):
        prev_id = ''
        cur_group, groups = [], []

        fns = sorted(os.listdir(self.source_fasta_dir))
        for fn in fns:
            cur_id = fn.split('_')[0]
            if cur_id != prev_id:
               groups.append(cur_group)
               cur_group = [fn.split('.')[0]]
               prev_id = cur_id
            else:
               cur_group.append(fn.split('.')[0])
        groups.append(cur_group)
        return groups[1:]

    def poly_g_link_all(self, fasta_groups):
        poly_g = 'G' * self.n_g
        chain_start_ids = {}
        for fasta_group in fasta_groups:
            utils.poly_g_link(self.source_fasta_dir, self.linked_fasta_dir, chain_start_ids, fasta_group, poly_g, self.n_g)
        return chain_start_ids

    def read_bd_resid_id_all(self):
        res = {}
        for id in self.pdb_ids:
            fn = join(self.input_pdb_dir, id + '.pdb')
            utils.read_chain_bd_resid_id_from_pdb(fn, id, res)
        return res

    def read_chain_names_all(self, gt_model_nm='native'):
        chain_names = {}
        for pdb_id in self.pdb_ids:
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            chain_name = utils.read_chain_name_from_pdb(fn)
            chain_names[pdb_id] = chain_name
        return chain_names

    def generate_fasta_from_pdb(self):
        ''' Generate source fasta based on gt pdb
              and polyg link all chains
            Also record residue id (1-based continuous between chains
              with polyg included) of start of each chain
        '''
        linker = 'G' * self.n_g
        valid_ids = []
        chain_start_resid_ids = {}

        for id in self.pdb_ids:
            in_fn = join(self.input_pdb_dir, id + '.pdb')
            seqs = utils.read_residue_from_pdb(in_fn)
            fasta = reduce(lambda acc, seq: acc + seq + linker, seqs, '')[:-self.n_g]
            if 'X' in fasta or 'x' in fasta:
                continue

            valid_ids.append(id)
            acc, start_ids = 1, []
            for seq in seqs:
                start_ids.append(acc)
                acc += len(seq) + self.n_g
            start_ids.append(len(fasta) + self.n_g + 1) # for ease of removal

            out_fn = join(self.linked_fasta_dir, id + '.fasta')
            utils.save_fasta(fasta, out_fn)
            chain_start_resid_ids[id] = start_ids

        np.save(self.pdb_ids_fn, np.array(valid_ids))
        return chain_start_resid_ids

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

            Input: chain_names   identifier for all chain of current gt pdb
                   chain_start_ids id of first residue in each chain (from pred pdb)
                   gt_chain_bd_ids id of two boundary residue of each chain (from gt pdb)
        '''
        if exists(pred_removed_linker_fn): return

        chain_names = self.chain_names[pdb_id]
        chain_start_ids = self.chain_start_ids[pdb_id]
        gt_chain_bd_ids = self.gt_chain_bd_ids[pdb_id]

        residues = utils.read_residue_from_pdb(pred_fn)[0]
        for i in range(1, len(chain_start_ids)-1):
            id = chain_start_ids[i] - 1 # 1-based to 0-based
            assert(residues[id - self.n_g : id] == 'G' * self.n_g)

        ppdb = PandasPdb()
        _ = ppdb.read_pdb(pred_fn)
        df = ppdb.df['ATOM']

        n = len(chain_start_ids)
        acc_linker_length = 0
        linker_atom_ids = []
        atom_los, atom_his, offsets = [], [], []

        for i in range(n-1): # for each chain
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
                           #(df['chain_id'] == chain_names[i]) ].index
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

            atom_los.append(atom_lo)
            atom_his.append(atom_hi)
            offsets.append(offset)
            df.loc[atom_lo:atom_hi ,'chain_id'] = chain_names[i]

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
        seqs1 = utils.read_residue_from_pdb(removed_linker_fn)
        seqs2 = utils.read_residue_from_pdb(gt_pdb_fn)
        ranges = utils.find_prune_ranges_all_chains(seqs1, seqs2, self.chain_names[pdb_id])
        utils.prune_renumber_seq_given_ranges(pred_removed_linker_fn, aligned_fn,
                                              self.chain_names[pdb_id], ranges[0], self.gt_chain_bd_ids[pdb_id])

    def calculate_rmsd(self, pdb_id, gt_pdb_fn, pred_pdb_fn, complex_fn, verbose=False):
        ''' Calculate rmsd between gt and pred
              superimpose pred onto gt (with receptor only)
              assume chain_names[-1] is the target chain (e.g. peptide or idr)
        '''
        rmsds = []
        utils.load_and_select\
            (self.interface_dist, gt_pdb_fn, pred_pdb_fn,
             self.chain_names[pdb_id], backbone=self.backbone,
             remove_hydrogen=self.remove_hydrogen)

        # superimpose receptor chains and calculate rmsd for idr
        super = cmd.super('native_R','pred_R')
        cmd.color('purple','native_R')
        cmd.color('yellow','native_L')
        cmd.color('gray','pred_R')
        cmd.color('orange','pred_L')
        cmd.multisave(complex_fn, 'all', format='pdb')

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

    def assert_fasta(args):
        ''' check if aa sequence from processed prediction pdb match sequence from gt pdb
        '''
        for id in args.pdb_ids:
            pdb_fn = join(args.input_pdb_dir, id + '.pdb')
            pred_fn = join(args.output_dir, id + '.fasta', args.pred_fn)
            seq_1 = utils.read_residue_from_pdb(pdb_fn)
            seq_2 = utils.read_residue_from_pdb(pred_fn)
            assert(seq_1 == seq_2)
