
import os
import json
import pickle
import numpy as np

from pymol import cmd
from functools import reduce
from os.path import join, exists

from utils import common as utils
from utils.dockq import calc_metrics


class pipeline:
    def __init__(self, option, args):
        self.verbose = args.verbose

        self.n_g = args.n_g
        self.input_pdb_dir = args.input_pdb_dir
        self.chain_ids_fn = args.chain_ids_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn

        if option == 'init_atom_prune':
            self._init_atom_prune(args)
        elif option == 'init_input_procs':
            self._init_input_processing(args)
        elif option == 'init_input_procs_fasta':
            self._init_input_processing_from_fasta(args)
        elif option == 'init_atom_locating':
            self._init_atom_locating(args)
        elif option == 'init_metric_plotting':
            self._init_metric_plotting(args)

    ########
    # inits
    def _init_input_processing(self, args):
        self.download_pdb = args.download_pdb
        self.pdb_ids_fn = args.pdb_ids_fn
        self.ds_spec_fn = args.dataset_spec_fn
        self.linked_fasta_dir = args.linked_fasta_dir
        self.pdb_download_command = args.pdb_download_command

    def _init_input_processing_from_fasta(self, args):
        self.download_pdb = args.download_pdb
        self.pdb_ids_fn = args.pdb_ids_fn
        self.ds_spec_fn = args.dataset_spec_fn
        self.source_fasta_dir = args.source_fasta_dir
        self.linked_fasta_dir = args.linked_fasta_dir
        self.pdb_download_command = args.pdb_download_command

    def _init_atom_prune(self, args):
        self.prune_pdb_ids = args.prune_pdb_ids
        self.atom_prune_ranges = args.atom_prune_ranges

    def _init_atom_locating(self, args):
        self.pdb_ids = args.pdb_ids

        self.from_fasta = args.from_fasta
        self.output_dir = args.output_dir
        self.pred_fn_str = args.pred_fn_str
        self.aligned_fn_str = args.aligned_fn_str

    #############
    # processors
    def process_input_base(self):
        self.parse_dataset_spec()
        if self.download_pdb:
            print("= Downloading source pdb")
            utils.run_bash(self.pdb_download_command)
        #print("= Downloading source fasta")
        #utils.run_bash(self.fasta_download_command)

        gt_chain_bd_resid_ids = self.read_bd_resid_id_all()
        with open(self.gt_chain_bd_resid_ids_fn, 'wb') as fp:
            pickle.dump(gt_chain_bd_resid_ids, fp)

        # get chain names and store locally
        #chain_ids = self.read_chain_ids_all()
        #with open(self.chain_ids_fn, 'wb') as fp:
        #    pickle.dump(chain_ids, fp)

    def process_input(self):
        self.process_input_base()
        chain_start_resid_ids = self.generate_fasta_from_pdb()
        with open(self.chain_start_resid_ids_fn, 'wb') as fp:
            pickle.dump(chain_start_resid_ids, fp)

    def process_input_from_fasta(self):
        self.process_input_base()
        groups = self.parse_to_groups()
        chain_start_resid_ids = self.poly_g_link_all(groups)
        with open(self.chain_start_resid_ids_fn, 'wb') as fp:
            pickle.dump(chain_start_resid_ids, fp)

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

    ##########
    # helpers
    def parse_dataset_spec(self):
        # read pdb_ids and correspd chain ids
        chain_ids, pdb_ids = {}, []

        f = open(self.ds_spec_fn)
        data = json.load(f)
        for pdb_id, v in data.items():
            pdb_ids.append(pdb_id)
            chain_ids[pdb_id] = [v['receptor_chain_ids'][0],
                                 v['idp_chain_ids'][0]]

        pdb_ids = np.array(pdb_ids)
        pdb_ids.sort()
        self.pdb_ids = pdb_ids
        pdb_id_str = reduce(lambda cur, acc: acc + ',' + cur, pdb_ids, '')

        # save pdb ids and chain ids
        # npy file for data processing
        np.save(self.pdb_ids_fn + '.npy', pdb_ids)
        # plain txt file for pdb downloading
        with open(self.pdb_ids_fn + '.txt', 'w') as fp:
            fp.write(pdb_id_str)
        with open(self.chain_ids_fn, 'wb') as fp:
            pickle.dump(chain_ids, fp)

    def parse_to_groups(self):
        '''
        '''
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

    def read_chain_ids_all(self, gt_model_nm='native'):
        chain_ids = {}
        for pdb_id in self.pdb_ids:
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            chain_id = utils.read_chain_id_from_pdb(fn)
            chain_id = utils.assign_receptor_ligand(fn, chain_id)
            chain_ids[pdb_id] = chain_id
        return chain_ids

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

    def assert_fasta(args):
        ''' check if aa sequence from processed prediction pdb match sequence from gt pdb
        '''
        for id in args.pdb_ids:
            pdb_fn = join(args.input_pdb_dir, id + '.pdb')
            pred_fn = join(args.output_dir, id + '.fasta', args.pred_fn)
            seq_1 = utils.read_residue_from_pdb(pdb_fn)
            seq_2 = utils.read_residue_from_pdb(pred_fn)
            assert(seq_1 == seq_2)
