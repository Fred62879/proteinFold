
import os
import json
import pickle
import numpy as np
import utils.common as utils
import utils.input_procs as iutils

from pymol import cmd
from functools import reduce
from os.path import join, exists
from utils.common import get_secondary_structure_resid_ids


class pipeline:
    def __init__(self, option, args):
        self.verbose = args.verbose

        self.n_g = args.n_g
        self.skip_fasta_download = args.skip_fasta_download
        self.generate_fasta_from_pdb = args.generate_fasta_from_pdb

        self.input_pdb_dir = args.input_pdb_dir
        self.orig_chain_ids_fn = args.orig_chain_ids_fn
        self.ordered_chain_ids_fn = args.ordered_chain_ids_fn
        self.chain_start_resid_ids_fn = args.chain_start_resid_ids_fn
        self.gt_chain_bd_resid_ids_fn = args.gt_chain_bd_resid_ids_fn
        self.second_structs_resid_ids_fn = args.second_structs_resid_ids_fn

        self.secondary_structures = args.secondary_structures

        if option == 'init_input_procs':
            self._init_input_processing(args)
        elif option == 'init_atom_locating':
            self._init_atom_locating(args)
        elif option == 'init_atom_prune':
            self._init_atom_prune(args)
        else:
            raise Exception('= !Unsupported input processing operation!')

    ########
    # inits
    def _init_input_processing(self, args):
        if self.generate_fasta_from_pdb:
            self.fasta_procs_func = self.extract_fasta_from_pdb
        else: self.fasta_procs_func = self.manually_download_fasta

        self.strategy = args.strategy
        self.download_pdb = args.download_pdb
        self.renumber = args.renumber_residues
        self.prune = args.prune_extra_residues
        self.has_dataset_spec = args.has_dataset_spec
        self.remove_linker = args.strategy == 'poly_g_link'

        self.pdb_ids_fn = args.pdb_ids_fn
        self.ds_spec_fn = args.dataset_spec_fn
        self.input_fasta_dir = args.input_fasta_dir
        self.source_fasta_dir = args.source_fasta_dir
        self.pdb_to_download_fn = args.pdb_to_download_fn
        self.pdb_download_command = args.pdb_download_command

    def _init_atom_prune(self, args):
        self.prune_pdb_ids = args.prune_pdb_ids
        self.atom_prune_ranges = args.atom_prune_ranges

    def _init_atom_locating(self, args):
        self.generate_fasta_from_pdb = args.generate_fasta_from_pdb
        self.output_dir = args.output_dir
        self.pred_fn_str = args.pred_fn_str
        self.aligned_fn_str = args.aligned_fn_str

    #############
    # processors
    def process_input(self):
        # If dataset spec file is provided, get pdb ids, chain ids from it
        # Otherwise, parse manually downloaded data to gather these
        if self.has_dataset_spec:
            self.parse_dataset_spec()
        else: self.parse_manual_data()

        # process fasta files
        if not self.skip_fasta_download or self.strategy == 'poly_g_link':
            self.fasta_procs_func()

        # use common set of pdb ids presented in pdb & fasta dir as experiment pdb ids
        self.collect_experiment_pdbs()

        # to remove linker or renumber residue id during output processing
        # we need id of boundary residue for each chain
        self.read_bd_resid_ids_all()

        # store original chain ids (order correspd. chain in gt pdb) & ordered chain ids
        self.read_chain_ids_all()

        # self.get_secondary_structure_ids()

    def locate_extra_atoms(self):
        for pdb_id in self.pdb_ids:
            print(f'{pdb_id} extra atoms: ')
            pdb_fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            dir = join(self.output_dir, pdb_id + '.fasta')
            if self.generate_fasta_from_pdb:
                pred_fn = join(dir, self.removed_linker_fn_str)
            else: pred_fn = join(dir, self.aligned_fn_str)
            ret = utils.find_residue_diff_in_atom_counts(pdb_fn, pred_fn)
            print(ret)

    def prune_extra_atoms(self):
        ''' prune_pdb_ids and prune_ranges are hardcoded in parser.py '''
        for pdb_id, rnge in zip(self.prune_pdb_ids, self.atom_prune_ranges):
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            iutils.prune_pdb_atoms(fn, rnge)

    ##########
    # helpers
    def parse_dataset_spec(self):
        ordered_chain_ids, pdb_ids = {}, []
        f = open(self.ds_spec_fn)
        data = json.load(f)
        for pdb_id, v in data.items():
            pdb_ids.append(pdb_id)
            ordered_chain_ids[pdb_id] = [v['receptor_chain_ids'][0],
                                         v['idp_chain_ids'][0]]
        with open(self.ordered_chain_ids_fn, 'wb') as fp:
            pickle.dump(ordered_chain_ids, fp)
        pdb_ids = np.array(pdb_ids)
        self.download_pdbs(pdb_ids)

        # for fasta downloading
        pdb_ids.sort()
        pdb_ids_str = reduce(lambda acc, cur: acc + ',' + cur, pdb_ids, '')
        with open(self.pdb_ids_fn + '.txt', 'w') as fp:
            fp.write(pdb_ids_str)

    def parse_manual_data(self):
        if self.download_pdb:
            fn = self.pdb_to_download_fn + '.txt'
            info = f'= Please provide a plain text file of comma separated pdb ids and store it here {fn}.\r\n'
            info += '  Press any key to continue when you finish.'
            input(info)
            assert(exists(fn))
            with open(fn) as fp:
                data = fp.readlines()
            if len(data) > 0:
                pdb_ids = np.array(data[0].split(','))
                pdb_ids[-1] = pdb_ids[-1][:4]
                self.download_pdbs(pdb_ids)
            else:
                pdb_ids = None

        else:
            pdb_ids = utils.parse_pdb_ids(self.input_pdb_dir, '.pdb')

        if pdb_ids is not None:
            # for fasta downloading
            pdb_ids.sort()
            pdb_ids_str = reduce(lambda acc, cur: acc + ',' + cur, pdb_ids, '')
            with open(self.pdb_ids_fn + '.txt', 'w') as fp:
                fp.write(pdb_ids_str)

    def manually_download_fasta(self):
        ''' Wait for confirmation from user that fasta files are stored in directed folder
            Provide a plain text file with comma separated pdb ids for ease of downloading
        '''
        fasta_fn = self.pdb_ids_fn + '.txt'
        info = f'= Please download fasta files to {self.source_fasta_dir}.\r\n'
        info += f'  You can use pdb ids provided in {fasta_fn} for downloading.\r\n'
        info += '  Press any key to continue when you finish.'
        input(info)
        assert(exists(self.source_fasta_dir) and
               len(os.listdir(self.source_fasta_dir)) != 0)
        if self.strategy == 'poly_g_link':
            iutils.poly_g_link_all(self.strategy, self.source_fasta_dir, self.input_fasta_dir,
                                  self.chain_start_resid_ids_fn, self.n_g)

    def extract_fasta_from_pdb(self):
        ''' Generate source fasta based on gt pdb.
            If poly glycine link strategy is used, also polyg link all
              chains and record id (1-based continuous between chains with polyg
              included) of 1st residue for each chain for later linker removal
        '''
        print('= generating fasta from pdb')

        if self.strategy == 'poly_g_link':
            linker = 'G' * self.n_g
            chain_start_resid_ids = {}

        for id in self.pdb_ids:
            in_fn = join(self.input_pdb_dir, id + '.pdb')
            seqs = utils.read_residue_from_pdb(in_fn)

            if self.strategy == 'poly_g_link':
                utils.poly_g_link_pdb(seqs, self.n_g, self.input_fasta_dir, chain_start_resid_ids)
            else:
                for i, seq in enumerate(seqs):
                    out_fn = join(self.input_fasta_dir, id + f'_{i}.fasta')
                    utils.save_fasta(seq, out_fn)

        if self.strategy == 'poly_g_link':
            with open(self.chain_start_resid_ids_fn, 'wb') as fp:
                pickle.dump(chain_start_resid_ids, fp)

    def collect_experiment_pdbs(self):
        # gather experiment pdb ids as common ids from pdb and fasta dir
        pdb_ids = utils.parse_pdb_ids(self.input_pdb_dir, '.pdb')
        if not self.skip_fasta_download:
            pdb_ids_1 = utils.parse_pdb_ids(self.input_fasta_dir, '.fasta')
            pdb_ids = np.array(list(set(pdb_ids).intersection(pdb_ids_1)))
        pdb_ids.sort()
        np.save(self.pdb_ids_fn + '.npy', pdb_ids)
        self.pdb_ids = pdb_ids

    def download_pdbs(self, pdb_ids):
        # download all specified pdbs. If pdb dir exists, only download pdb not presented in the dir
        pdb_downloaded = utils.parse_pdb_ids(self.input_pdb_dir, '.pdb')
        pdb_to_download = np.array(list(set(pdb_ids) - set(pdb_downloaded)))
        print(pdb_downloaded, pdb_ids, pdb_to_download)
        pdb_to_download.sort()

        # save pdb ids as a str in text file for pdb downloading
        pdb_to_download_str = reduce(lambda acc, cur: acc + ',' + cur, pdb_to_download, '')
        with open(self.pdb_to_download_fn + '.txt', 'w') as fp:
            fp.write(pdb_to_download_str)

        print("= Downloading source pdb")
        print(f'- {len(pdb_downloaded)} pdbs already downloaded, {len(pdb_to_download)} pdbs to download')
        utils.run_bash(self.pdb_download_command)

    def read_bd_resid_ids_all(self):
        # read id of boundary (1st and last) residues of each chain
        bypass = (not self.remove_linker and not self.renumber and not self.prune) \
            or exists(self.gt_chain_bd_resid_ids_fn)
        if bypass: return

        print('= loading boundary residue ids')
        gt_chain_bd_resid_ids = {}
        for id in self.pdb_ids:
            fn = join(self.input_pdb_dir, id + '.pdb')
            utils.read_chain_bd_resid_id_from_pdb(fn, id, gt_chain_bd_resid_ids)
        with open(self.gt_chain_bd_resid_ids_fn, 'wb') as fp:
            pickle.dump(gt_chain_bd_resid_ids, fp)

    def read_chain_ids_all(self):
        ''' Read chain id (orig order and receptor-ligand order) for each pdb
            If we have dataset spec file, ordered chain ids are already obtained
            Otherwise, ordered chain ids are obtained via comparing num of atoms
              in each chain. Longer chain is receptor (1st), the other is ligand (2nd)
        '''
        bypass = exists(self.orig_chain_ids_fn) and \
            exists(self.ordered_chain_ids_fn)
        if bypass: return

        print('= loading chain ids')
        orig_chain_ids, ordered_chain_ids = {}, {}

        for pdb_id in self.pdb_ids:
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            orig_chain_id = utils.read_chain_id_from_pdb(fn)
            orig_chain_ids[pdb_id] = orig_chain_id
            if not self.has_dataset_spec:
                ordered_chain_id = utils.assign_receptor_ligand(fn, orig_chain_id)
                ordered_chain_ids[pdb_id] = ordered_chain_id

        with open(self.orig_chain_ids_fn, 'wb') as fp:
            pickle.dump(orig_chain_ids, fp)
        if not self.has_dataset_spec:
            with open(self.ordered_chain_ids_fn, 'wb') as fp:
                pickle.dump(ordered_chain_ids, fp)

    def get_secondary_structure_ids(self):
        second_structs_resid_ids = {}
        for pdb_id in self.pdb_ids:
            fn = join(self.input_pdb_dir, pdb_id + '.pdb')
            cur_resid_id = get_secondary_structure_resid_ids(fn, self.secondary_structures)
            second_structs_resid_ids[pdb_id] = cur_resid_id
        with open(self.second_structs_resid_ids_fn, 'wb') as fp:
            pickle.dump(second_structs_resid_ids, fp)
