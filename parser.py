import os
import numpy as np

from pathlib import Path
from os.path import join, exists


def add_cmd_line_args(parser):
    parser.add('-c', '--config', required=False, is_config_file=True)

    # experiment setup
    parser.add_argument('--strategy', type=str, help='poly_g_link or multimer')
    parser.add_argument('--pdb_cho', type=str, default='0', help='can ignore this')
    parser.add_argument('--operations', type=str, nargs='+', help='processing options, refer to config file')

    parser.add_argument('--code_dir', type=str, required=True, help='code directory')
    parser.add_argument('--data_dir', type=str, required=True, help='data directory')

    parser.add_argument('--model_name', type=str, default='alphafold',
                        help='differentiate between different configuration of alphafold, define as whatever you want')
    parser.add_argument('--dataset_name', type=str, default='ds0',
                        help='differentiate between datasets, define as whatever you want')
    parser.add_argument('--dataset_spec_name', type=str, default='rg-dimers',
                        help='name of json file containing input dataset information')

    parser.add_argument('--n_seqs', type=int, default=2, help='number of sequences to predict')
    parser.add_argument('--n_g', type=int, default=6, help='number of glycines thta compose the linker')
    parser.add_argument('--interface_dist', type=int, default=6, help='dist in angstrom to define interface residues')

    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--generate_fasta_from_pdb', action='store_true', default=False,
                        help='extra sequence from pdb as fasta or download fasta directly from rcsb.org')

    parser.add_argument('--download_pdb', action='store_true', default=False,
                        help='manually download all pdb files as specified by dataset spec file')
    parser.add_argument('--prune_pdb_atoms', action='store_true', default=False,
                        help='remove extra atoms from gt pdb file')
    parser.add_argument('--has_dataset_spec', action='store_true', default=False,
                        help='whether or not dataset spec file is provided')
    parser.add_argument('--skip_fasta_download', action='store_true', default=False,
                        help='skip fasta download step, can be True if multimer strategy is used')

    parser.add_argument('--dockq', action='store_true', default=False,
                        help='use official dockq code for rmsd and dockq calculation')
    parser.add_argument('--backbone', action='store_true', default=False,
                        help='calculate and visualize using only backbone residues')
    parser.add_argument('--remove_hydrogen', action='store_true', default=False)
    parser.add_argument('--renumber_residues', action='store_true', default=False,
                        help='renumber residue ids to be the same as gt pdb')
    parser.add_argument('--prune_extra_residues', action='store_true', default=False,
                        help='prune extra residues from alphafold prediction')


    args = parser.parse_args()
    config = vars(args)
    return config


def add_hardcoded_args(config):
    config['gt_model_nm'] = 'native'

    config['pdb_str'] = 'pdbs'
    config['input_str'] = 'input'
    config['output_str'] = 'output'

    config['source_fasta_dir_str'] = 'source_fasta'

    if config['strategy'] == 'poly_g_link':
        config['input_fasta_dir_str'] = 'poly_g_' + str(config['n_g'])
        config['output_dir_str'] = config['input_fasta_dir_str']
    else: config['output_dir_str'] = config['strategy']

    if config['generate_fasta_from_pdb']:
        config['input_fasta_dir_str'] += '_from_pdb'

    config['fnat_path_str'] = 'utils'
    config['pred_pdb_fn_str'] = 'ranked_0.pdb'
    config['ranking_fn_str'] = 'ranking_debug.json'
    config['pruned_fn_str'] = 'ranked_0_pruned.pdb'
    config['aligned_fn_str'] = 'ranked_0_aligned.pdb'
    config['removed_linker_fn_str'] = 'ranked_0_removed_linker.pdb'

    config['pdb_ids_fn_str'] = 'pdb_ids' # + config['pdb_cho']
    config['failed_pdbs_fn_str'] = 'failed_pdbs.pkl'
    config['pdb_exclude_fn_str'] = 'pdb_exclude' #+ config['pdb_cho']
    config['pdb_gpu_done_fn_str'] = 'pdb_gpu_done' #+ config['pdb_cho']
    config['pdb_to_download_fn_str'] = 'pdb_to_download' # + config['pdb_cho']
    config['orig_chain_ids_fn_str'] = 'orig_chain_ids.pkl' #+config['pdb_cho']+'.pkl'
    config['ordered_chain_ids_fn_str'] = 'ordered_chain_ids.pkl' #+config['pdb_cho']+'.pkl'
    config['chain_start_resid_ids_fn_str'] = 'chain_start_ids.pkl' #+config['pdb_cho']+'.pkl'
    config['gt_chain_bd_resid_ids_fn_str'] = 'gt_chain_bd_ids.pkl' #+config['pdb_cho']+'.pkl'

    backbone_str = '' if not config['backbone'] else '_backbone'
    config['metric_fn_str'] = f'metric{backbone_str}'
    if config['dockq']:
        config['metric_col_names'] = ['pdb_id','irms','Lrms','dockQ']
        config['metric_fn_str'] += '_dockq'
    else:
        config['metric_col_names'] = ['pdb_id','ligand (super r)','interface ligand (super r)']
    config['metric_col_names'].extend(['idr_len','num_interface_resid','plddt','sasa'])
    config['metric_fn_str'] += '.csv'

    # add range of atoms to prune
    # hardcoded, corresponding to pdb_ids
    config['prune_pdb_ids'] = ['3CAA']

    # range is upper exclusive, 0-based
    # atom id in pdb file is 1-based
    # so need to subtract one
    config['atom_prune_ranges'] = [ [[[1602,1608]]] ]

def add_path(config):
    input_dir = join(config['data_dir'], config['input_str'])
    input_ds_dir = join(input_dir, config['dataset_name'])
    config['input_pdb_dir'] = join(input_ds_dir, config['pdb_str'])
    config['source_fasta_dir'] = join(input_ds_dir, config['source_fasta_dir_str'])
    config['output_dir'] = join(config['data_dir'], config['output_str'], config['dataset_name'] +
                                '_' + config['model_name'], config['output_dir_str'])

    if config['strategy'] == 'multimer':
        config['input_fasta_dir'] = config['source_fasta_dir']
    else: config['input_fasta_dir'] = join(input_ds_dir, config['input_fasta_dir_str'])

    # create dir
    for dir in [input_dir, input_ds_dir, config['output_dir'], config['input_pdb_dir'],
                config['input_fasta_dir']]:
        if not exists(dir):
            Path(dir).mkdir(parents=True, exist_ok=True)

    # add filenames
    config['pdb_ids_fn'] = join(input_ds_dir, config['pdb_ids_fn_str'])
    config['fnat_path'] = join(config['code_dir'], config['fnat_path_str'])
    config['dataset_spec_fn'] = join(input_dir, config['dataset_spec_name'])
    config['metric_fn'] = join(config['output_dir'], config['metric_fn_str'])
    config['dataset_info_fn'] = join(input_ds_dir, config['dataset_spec_name'])
    config['pdb_exclude_fn'] = join(input_ds_dir, config['pdb_exclude_fn_str']+'.npy')
    config['failed_pdbs_fn'] = join(config['output_dir'], config['failed_pdbs_fn_str'])
    config['pdb_to_download_fn'] = join(input_ds_dir, config['pdb_to_download_fn_str'])
    config['pdb_gpu_done_fn'] = join(input_ds_dir, config['pdb_gpu_done_fn_str']+'.npy')
    config['orig_chain_ids_fn'] = join(input_ds_dir, config['orig_chain_ids_fn_str'])
    config['ordered_chain_ids_fn'] = join(input_ds_dir, config['ordered_chain_ids_fn_str'])
    config['chain_start_resid_ids_fn'] = join(config['input_fasta_dir'], config['chain_start_resid_ids_fn_str'])
    config['gt_chain_bd_resid_ids_fn'] = join(config['input_fasta_dir'], config['gt_chain_bd_resid_ids_fn_str'])

def add_bash_commands(config):
    # bash command for data downloading
    download_script_fn = config['code_dir'] + "/scripts/batch_download.sh"
    pdb_to_download_fn = config['pdb_to_download_fn'] + '.txt'
    pdb_dir = config['input_pdb_dir']
    config['pdb_download_command'] = f"{download_script_fn} -f {pdb_to_download_fn} -o {pdb_dir} -p"

def assert_args(config):
    assert(not config['skip_fasta_download'] or config['strategy'] == 'multimer')

def parse_args(parser):
    print('=== Parsing ===')
    config = add_cmd_line_args(parser)
    add_hardcoded_args(config)
    add_path(config)
    add_bash_commands(config)
    return config
