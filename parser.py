import os
import numpy as np

from pathlib import Path
from os.path import join, exists


def add_cmd_line_args(parser):
    parser.add('-c', '--config', required=False, is_config_file=True)

    # experiment setup
    parser.add_argument('--pdb_cho', type=str, default='0')
    parser.add_argument('--operations', type=str, nargs='+')

    parser.add_argument('--code_dir', type=str, required=True)
    parser.add_argument('--data_dir', type=str, required=True)

    parser.add_argument('--model_name', type=str, default='alphafold')
    parser.add_argument('--dataset_name', type=str, default='ds0')
    parser.add_argument('--dataset_spec_name', type=str, default='rg-dimers')

    parser.add_argument('--n_seqs', type=int, default=2, help='number of sequences to predict')
    parser.add_argument('--n_g', type=int, default=6, help='number of glycines thta compose the linker')
    parser.add_argument('--interface_dist', type=int, default=6, help='dist in angstrom to define interface residues')

    parser.add_argument('--dockq', action='store_true', default=False)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--backbone', action='store_true', default=False)
    parser.add_argument('--from_fasta', action='store_true', default=False)
    parser.add_argument('--download_pdb', action='store_true', default=False)
    parser.add_argument('--prune_pdb_atoms', action='store_true', default=False,
                        help='remove extra atoms from pdb file, currently only used for native pdb')
    parser.add_argument('--prune_and_renumber', action='store_true', default=False)
    parser.add_argument('--separate_fasta', action='store_true', default=False,
                        help='whether each monomer is stored in a single file or all monomers store in same fasta file')
    parser.add_argument('--prune_unknown', action='store_true', default=False)
    parser.add_argument('--remove_hydrogen', action='store_true', default=False)

    args = parser.parse_args()
    config = vars(args)
    return config


def add_hardcoded_args(config):
    #config['code_dir'] = config['code_dir'][1:-1]
    #config['data_dir'] = config['data_dir'][1:-1]

    config['gt_model_nm'] = 'native'

    config['pdb_str'] = 'pdbs'
    config['input_str'] = 'input'
    config['output_str'] = 'output'

    config['source_fasta_dir_str'] = 'source_fasta'
    config['linked_seq_dir_str'] = 'poly_g_' + str(config['n_g'])
    if config['from_fasta']:
        config['linked_seq_dir_str'] += '_fasta'

    config['fnat_path_str'] = 'utils'
    config['pred_fn_str'] = 'ranked_0.pdb'
    config['complex_fn_str'] = 'complex.pdb'
    config['ranking_fn_str'] = 'ranking_debug.json'
    config['removed_linker_fn_str'] = 'ranked_0_removed_linker.pdb'
    config['aligned_fn_str'] = 'ranked_0_removed_linker_aligned.pdb'

    config['pdb_ids_fn_str'] = 'pdb_ids' # + config['pdb_cho']
    config['pdb_gpu_done_fn_str'] = 'pdb_gpu_done' #+ config['pdb_cho']
    config['pdb_exclude_fn_str'] = 'pdb_exclude' #+ config['pdb_cho']
    config['chain_ids_fn_str'] = 'chain_ids.pkl' #+config['pdb_cho']+'.pkl'
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
    config['linked_fasta_dir'] = join(input_ds_dir, config['linked_seq_dir_str'])
    config['source_fasta_dir'] = join(input_ds_dir, config['source_fasta_dir_str'])
    config['output_dir'] = join(config['data_dir'], config['output_str'], config['dataset_name'] +
                                '_' + config['model_name'], config['linked_seq_dir_str'])
    # create dir
    for dir in [input_dir, input_ds_dir, config['output_dir'], config['input_pdb_dir'],
                config['linked_fasta_dir'], config['linked_fasta_dir']]:
        if not exists(dir):
            Path(dir).mkdir(parents=True, exist_ok=True)

    # add filenames
    config['pdb_ids_fn'] = join(input_ds_dir, config['pdb_ids_fn_str'])
    config['fnat_path'] = join(config['code_dir'], config['fnat_path_str'])
    config['dataset_spec_fn'] = join(input_dir, config['dataset_spec_name'])
    config['dataset_info_fn'] = join(input_ds_dir, config['dataset_spec_name'])
    config['pdb_gpu_done_fn'] = join(input_ds_dir, config['pdb_gpu_done_fn_str']+'.npy')
    config['pdb_exclude_fn'] = join(input_ds_dir, config['pdb_exclude_fn_str']+'.npy')
    config['metric_fn'] = join(config['output_dir'], config['metric_fn_str'])
    config['chain_ids_fn'] = join(config['linked_fasta_dir'], config['chain_ids_fn_str'])
    config['chain_start_resid_ids_fn'] = join(config['linked_fasta_dir'], config['chain_start_resid_ids_fn_str'])
    config['gt_chain_bd_resid_ids_fn'] = join(config['linked_fasta_dir'], config['gt_chain_bd_resid_ids_fn_str'])

def add_bash_commands(config):
    # bash command for data downloading
    download_script_fn = config['code_dir'] + "/scripts/batch_download.sh"
    pdb_download_fn = config['pdb_ids_fn'] + '.txt'
    pdb_dir = config['input_pdb_dir']
    config['pdb_download_command'] = f"{download_script_fn} -f {pdb_download_fn} -o {pdb_dir} -p"

def parse_args(parser):
    print('=== Parsing ===')
    config = add_cmd_line_args(parser)
    add_hardcoded_args(config)
    add_path(config)
    add_bash_commands(config)
    return config
