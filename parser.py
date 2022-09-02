import os
import numpy as np

from os.path import join, exists


def add_cmd_line_args(parser):
    parser.add('-c', '--config', required=False, is_config_file=True)

    # experiment setup
    parser.add_argument('--operations', type=str, nargs='+')

    parser.add_argument('--trail_id', type=str, default='trail_dum')
    parser.add_argument('--model_name', type=str, default='model_dum')
    parser.add_argument('--experiment_id', type=str, default='exp_dum')

    parser.add_argument('--pred_fn', type=str, default='ranked_0')

    parser.add_argument('--remove_x', action='store_true', default=False)
    parser.add_argument('--backbone', action='store_true', default=False)
    parser.add_argument('--prune_unknown', action='store_true', default=False)
    parser.add_argument('--remove_hydrogen', action='store_true', default=False)

    args = parser.parse_args()
    config = vars(args)
    return config

def add_hardcoded_args(config):
    config['n_g'] = 6
    config['gt_model_nm'] = 'native'

    config['data_dir'] = '/media/fred/Local Disk/Projects/bioinfo/data'

    config['pdb_str'] = 'pdbs'
    config['input_str'] = 'input'
    config['output_str'] = 'output'
    config['source_dir_str'] = 'source'
    config['linker_dir_str'] = 'poly_g_6'

    config['reorder_fn_str'] = 'reordered'
    config['pred_fn_str'] = 'ranked_0.pdb'
    config['complex_fn_str'] = 'complex.pdb'
    config['chain_names_fn_str'] = 'chain_names.pkl'
    backbone_str = '' if not config['backbone'] else '_backbone'
    config['rmsd_fn_str'] = f'rmsd{backbone_str}.csv'
    config['chain_start_resid_ids_fn_str'] = 'chain_start_ids.pkl'
    config['gt_chain_bd_resid_ids_fn_str'] = 'gt_chain_bd_ids.pkl'
    config['removed_linker_fn_str'] = 'ranked_0_removed_linker.pdb'

    #config['pdb_ids'] = ['2QTV']
    config['pdb_ids'] = ['1YCQ','2AZE','2RSN']
    #config['pdb_ids'] = ['1AWR','1EG4','1ELW','1ER8','1JD5']
    #config['pdb_ids'] = ['1YCQ','2AZE','2M3M','2QTV','2RSN','3DF0','4U7T']
    config['rmsd_names'] = ['pdb_id','rmsd_init','rmsd_idr_super_r','rmsd_r_super_r','rmsd_r_super_idr','rmsd_idr_align','rmsd_r_align']


def add_path(config):
    input_dir = join(config['data_dir'], config['input_str'], config['experiment_id'])

    config['input_pdb_dir'] = join(input_dir, config['pdb_str'])
    config['source_fasta_dir'] = join(input_dir, config['source_dir_str'])
    config['linker_fasta_dir'] = join(input_dir, config['linker_dir_str'])
    config['output_dir'] = join(config['data_dir'], config['output_str'], config['experiment_id'] + '_' + config['model_name'], config['linker_dir_str'])

    config['rmsd_fn'] = join(config['output_dir'], config['rmsd_fn_str'])
    config['chain_names_fn'] = join(config['linker_fasta_dir'], config['chain_names_fn_str'])
    config['chain_start_resid_ids_fn'] = join(config['linker_fasta_dir'], config['chain_start_resid_ids_fn_str'])
    config['gt_chain_bd_resid_ids_fn'] = join(config['linker_fasta_dir'], config['gt_chain_bd_resid_ids_fn_str'])


def parse_args(parser):
    print('=== Parsing ===')
    config = add_cmd_line_args(parser)
    add_hardcoded_args(config)
    add_path(config)
    return config
