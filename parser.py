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
    parser.add_argument('--trail_id', type=str, default='trail_dum')
    parser.add_argument('--model_name', type=str, default='model_dum')
    parser.add_argument('--experiment_id', type=str, default='exp_dum')
    parser.add_argument('--n_g', type=int, default=6, help='number of glycines thta compose the linker')
    parser.add_argument('--n_seqs', type=int, default=2, help='number of sequences to predict')
    parser.add_argument('--interface_dist', type=int, default=6, help='dist in angstrom to define interface residues')

    parser.add_argument('--pred_fn', type=str, default='ranked_0')

    parser.add_argument('--dockq', action='store_true', default=False)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--remove_x', action='store_true', default=False)
    parser.add_argument('--backbone', action='store_true', default=False)
    parser.add_argument('--from_fasta', action='store_true', default=False)
    parser.add_argument('--prune_and_renumber', action='store_true', default=False)
    parser.add_argument('--separate_fasta', action='store_true', default=False,
                        help='whether each monomer is stored in a single file or all monomers store in same fasta file')
    parser.add_argument('--prune_unknown', action='store_true', default=False)
    parser.add_argument('--remove_hydrogen', action='store_true', default=False)

    args = parser.parse_args()
    config = vars(args)
    return config


def add_hardcoded_args(config):
    config['gt_model_nm'] = 'native'

    config['pdb_str'] = 'pdbs'
    config['input_str'] = 'input'
    config['output_str'] = 'output'

    config['source_fasta_dir_str'] = 'source_fasta'
    config['linked_seq_dir_str'] = 'poly_g_' + str(config['n_g'])
    if config['from_fasta']:
        config['linked_seq_dir_str'] += '_fasta'

    config['pred_fn_str'] = 'ranked_0.pdb'
    config['complex_fn_str'] = 'complex.pdb'
    config['removed_linker_fn_str'] = 'ranked_0_removed_linker.pdb'
    config['aligned_fn_str'] = 'ranked_0_removed_linker_aligned.pdb'

    config['pdb_ids_fn_str'] = 'pdb_ids' + config['pdb_cho']
    config['pdb_gpu_done_fn_str'] = 'pdb_gpu_done' + config['pdb_cho']
    config['pdb_exclude_fn_str'] = 'pdb_exclude' + config['pdb_cho']
    config['chain_names_fn_str'] = 'chain_names'+config['pdb_cho']+'.pkl'
    config['chain_start_resid_ids_fn_str'] = 'chain_start_ids'+config['pdb_cho']+'.pkl'
    config['gt_chain_bd_resid_ids_fn_str'] = 'gt_chain_bd_ids'+config['pdb_cho']+'.pkl'

    backbone_str = '' if not config['backbone'] else '_backbone'
    config['rmsd_fn_str'] = f'rmsd{backbone_str}'
    if config['dockq']:
        config['rmsd_names'] = ['pdb_id','irms','Lrms','dockQ']
        config['rmsd_fn_str'] += '_dockq'
    else:
        config['rmsd_names'] = ['pdb_id','ligand (super r)','interface ligand (super r)']
    config['rmsd_fn_str'] += '.csv'

def add_path(config):
    input_dir = join(config['data_dir'], config['input_str'], config['experiment_id'])

    config['input_pdb_dir'] = join(input_dir, config['pdb_str'])
    config['linked_fasta_dir'] = join(input_dir, config['linked_seq_dir_str'])
    config['source_fasta_dir'] = join(input_dir, config['source_fasta_dir_str'])
    config['output_dir'] = join(config['data_dir'], config['output_str'], config['experiment_id'] + '_' + config['model_name'], config['linked_seq_dir_str'])

    # create dir
    for dir in [input_dir, config['output_dir'], config['input_pdb_dir'],
                config['linked_fasta_dir'], config['linked_fasta_dir']]:
        if not exists(dir):
            Path(dir).mkdir(parents=True, exist_ok=True)

    # add filenames
    config['pdb_ids_fn'] = join(input_dir, config['pdb_ids_fn_str']+'.npy')
    config['pdb_gpu_done_fn'] = join(input_dir, config['pdb_gpu_done_fn_str']+'.npy')
    config['pdb_exclude_fn'] = join(input_dir, config['pdb_exclude_fn_str']+'.npy')
    config['rmsd_fn'] = join(config['output_dir'], config['rmsd_fn_str'])
    config['chain_names_fn'] = join(config['linked_fasta_dir'], config['chain_names_fn_str'])
    config['chain_start_resid_ids_fn'] = join(config['linked_fasta_dir'], config['chain_start_resid_ids_fn_str'])
    config['gt_chain_bd_resid_ids_fn'] = join(config['linked_fasta_dir'], config['gt_chain_bd_resid_ids_fn_str'])


def select_pdb_ids(config):
    #config['pdb_ids'] = ['2XZE']
    #config['pdb_ids'] = ['1YCQ','2AZE'] #,'2RSN']
    #config['pdb_ids'] = ['1AWR','1EG4','1ELW','1ER8','1JD5']
    #config['pdb_ids'] = ['1YCQ','2AZE','2M3M','2QTV','2RSN','3DF0','4U7T']

    #fn = config['pdb_ids_fn']
    fn = config['pdb_gpu_done_fn']
    excld_fn = config['pdb_exclude_fn']

    if exists(fn):
        pdb_ids = np.load(fn)
        if exists(excld_fn):
            exclude_ids = np.load(excld_fn)
            #exclude_ids = np.append(exclude_ids,['2M3M','4X34'])
            pdb_ids = np.array(list( set(pdb_ids) - set(exclude_ids) ))
        pdb_ids = np.sort(pdb_ids)
    else:
        #if config['from_fasta']:
        #    all_pdb_ids = np.array(os.listdir(config['source_fasta_dir']))
        #else:
        all_pdb_ids = np.array(os.listdir(config['input_pdb_dir']))
        pdb_ids = np.random.choice(all_pdb_ids, config['n_seqs'], replace=False)
        split_char = '_' if config['from_fasta'] and config['separate_fasta'] else '.'
        pdb_ids = np.array([id.split(split_char)[0] for id in pdb_ids])
        pdb_ids = np.sort(pdb_ids)
        np.save(fn, pdb_ids)

    config['pdb_ids'] = pdb_ids[:1]
    print(f'selected proteins {pdb_ids}')


def parse_args(parser):
    print('=== Parsing ===')
    config = add_cmd_line_args(parser)
    add_hardcoded_args(config)
    add_path(config)
    select_pdb_ids(config)
    return config
