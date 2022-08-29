
import os
import pandas
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from pymol import cmd
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
from functools import reduce
from biopandas.pdb import PandasPdb
from Bio.SeqRecord import SeqRecord


def read_fasta(fn):
    fasta_sequences = SeqIO.parse(open(fn),'fasta')
    for fasta in fasta_sequences:
        monomer = str(fasta.seq)
        break
    return monomer

def read_residue_from_pdb(fn, print=False):
    fastas = []
    with open(fn, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            if print:
                print('>' + record.id)
                print(record.seq)
            else:
                _ = record.id
                fastas.append(str(record.seq))
    return fastas

def generate_fasta_from_pdb(in_dir, out_dir, pdb_ids, n_g):
    ''' generate source fasta based on gt pdb
          and polyg link all chains
        also record residue id (1-based continuous between chains
          with polyg included) of start of each chain
    '''
    linker = 'G'*n_g
    chain_start_resid_ids = {}

    for id in pdb_ids:
        in_fn = join(in_dir, id + '.pdb')
        seqs = read_residue_from_pdb(in_fn)
        fasta = reduce(lambda acc, seq: acc + seq + linker, seqs, '')[:-n_g]

        acc, start_ids = 1, []
        for seq in seqs:
            start_ids.append(acc)
            acc += len(seq) + n_g
        start_ids.append(len(fasta) + n_g + 1) # for ease of removal

        out_fn = join(out_dir, id + '.fasta')
        save_fasta(fasta, out_fn)
        chain_start_resid_ids[id] = start_ids
        print(start_ids)
    return chain_start_resid_ids

def save_fasta(seq, fn, id='0'):
    seq = SeqRecord(Seq(seq),id=id)
    SeqIO.write(seq, fn, 'fasta')

def read_chain_bd_resid_id_from_pdb(pdb_fn, pdb_id, bd_resid_ids):
    ''' return a list containing id of two boundary
          residues in each chain
    '''

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    cur_bd_resid_ids = []
    chain_ids = read_chain_id_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        cur_bd_resid_ids.append([ dup_resid_ids[0], dup_resid_ids[-1] ])

    bd_resid_ids[pdb_id] = cur_bd_resid_ids
    return bd_resid_ids

def read_bd_resid_id_all(dir, pdb_ids):
    res = {}
    for id in pdb_ids:
        fn = join(dir, id + '.pdb')
        read_chain_bd_resid_id_from_pdb(fn, id, res)
    return res

def read_residue_id_from_pdb(pdb_fn):
    ''' return a list containing unique residue ids for the whole pdb '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']

    resid_ids = []
    chain_ids = read_chain_id_from_pdb(pdb_fn)
    for chain_id in chain_ids:
        cur_chain_atom_ids = df['chain_id'] == chain_id
        dup_resid_ids = list(df.loc[cur_chain_atom_ids,'residue_number'])
        lo, hi = dup_resid_ids[0], dup_resid_ids[-1]
        de_dup_resid_ids = list(np.arange(lo,hi+1))
        resid_ids.extend(de_dup_resid_ids)

    return resid_ids

def read_chain_id_from_pdb(pdb_fn):
    ''' read chain identifier, de-depulicate while preserving order '''
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']
    return list(dict.fromkeys(df.loc[:,'chain_id']))

def count_num_atoms_each_residue(pdb_fn):
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']
    n = df.shape[0]

    num_atoms = []
    start_atom_id = 0
    prev_id = df.loc[0,'residue_number']
    resid_ids = [prev_id]

    for atom_id in range(1, n):
       cur_id = df.loc[atom_id, 'residue_number']
       if cur_id != prev_id:
           num_atoms.append(atom_id - start_atom_id)
           start_atom_id = atom_id
           prev_id = cur_id

    num_atoms.append(n - start_atom_id)
    return num_atoms

def get_chain_identifiers_all(pdb_ids, input_pdb_dir, gt_model_nm='native'):
    chain_names = {}
    for pdb_id in pdb_ids:
        fn = join(input_pdb_dir, pdb_id + '.pdb')
        #chain_name = split_chains_wo_env(fn, model=gt_model_nm)
        chain_name = read_chain_id_from_pdb(fn)
        chain_names[pdb_id] = chain_name
    return chain_names

def remove_linker(pred_fn, out_fn, linker_len, chain_names,
                  chain_start_ids, gt_chain_bd_ids):

    ''' remove linker from predicted pdb
          also reorder residue id and rename chain identifier
          as per gt pdb

        Note: pred pdb residue id is continuous 1-indexed
          and since all chains are polyg linked, resid id
          is continuous between chains (polyg linker inclusive)

        chain_names:     identifier for all chain of current gt pdb
        chain_start_ids: id of first residue in each chain (from pred pdb)
        gt_chain_bd_ids: id of two boundary residue of each chain (from gt pdb)
    '''

    residues = read_residue_from_pdb(pred_fn)[0]
    for i in range(1, len(chain_start_ids)-1):
        id = chain_start_ids[i] - 1 # 1-based to 0-based
        print(residues[id-linker_len:id])
        assert(residues[id-linker_len : id] == 'GGGGGG')

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pred_fn)
    df = ppdb.df['ATOM']

    n = len(chain_start_ids)
    linker_atom_ids = []
    atom_los, atom_his, offsets = [], [], []

    for i in range(n-1): # for each chain
        resid_lo = chain_start_ids[i]
        # resid_hi is residue id of first 'g' in pred
        resid_hi = chain_start_ids[i+1] - linker_len

        # make sure gt and pred chain has same length
        gt_resid_lo, gt_resid_hi = gt_chain_bd_ids[i]
        assert(gt_resid_hi - gt_resid_lo + 1 == resid_hi - resid_lo)

        # locate all atoms in pred pdb for current chain
        atom_ids = df[ (df['residue_number'] >= resid_lo) &
                       (df['residue_number'] < resid_hi) ].index
        atom_lo, atom_hi = atom_ids[0], atom_ids[-1]

        ''' reorder residue id and rename chain
            df.loc upper bound is inclusive
            atom_hi is id of last atom of current chain
        '''
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

    # renumber chain residue id
    for atom_lo, atom_hi, offset in zip(atom_los, atom_his, offsets):
        df.loc[atom_lo:atom_hi ,'residue_number'] -= offset
        print(df.loc[atom_lo,'residue_number'],df.loc[atom_hi,'residue_number'])

    # drop linker residues
    df.drop(linker_atom_ids, axis=0, inplace=True)
    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def calculate_rmsd(gt_pdb_fn, pred_pdb_fn, chain_names):
    ''' calculate rmsd between gt and pred

        superimpose pred onto gt (with receptor only)
          assume chain_names[0] is the receptor chain name
    '''
    cmd.load(gt_pdb_fn, 'native')
    cmd.load(pred_pdb_fn, 'pred')

    r_name, o_name = chain_names
    cmd.select('pred_R', f'pred and chain {r_name}')
    cmd.select('pred_O', f'pred and chain {o_name}')
    cmd.select('native_R', f'native and chain {r_name}')
    cmd.select('native_O', f'native and chain {o_name}')

    metrics = cmd.align('native_O','pred_O')
    print(metrics)

    # superimpose pred onto gt
    super = cmd.super('pred_R','native_R')
    metrics = cmd.align('pred_O','native_O')
    print(metrics)


def reorder_csv(in_file, out_file, colnames):
    with open(in_file, 'r') as infile, open(out_file, 'a') as outfile:
        # output dict needs a list for new column ordering
        writer = csv.DictWriter(outfile, fieldnames=colnames)
        # reorder the header first
        writer.writeheader()
        for row in csv.DictReader(infile):
            # writes the reordered rows to the new file
            writer.writerow(row)

def peptide_metric(pdb_ids, input_dir, output_dir, pred_fn, chain_names):
    metrics = []

    for id in pdb_ids[0:1]:
    #for model_file in glob.glob(os.path.join(input_dir, 'linker_removed*.pdb')):
	#model_name = os.path.basename(model_file).replace('.pdb', '')
	#rank = int(str(rank_re.search(model_name).group(0)).replace('rank_', ''))
	#model_no = int(str(model_re.search(model_name).group(0)).replace('model_', ''))
        model_name = 'afold'
        rank = 0
        model_no = 'dum'
        rec_chain, pep_chain = chain_names[id]
        metrics_for_a_model = [model_name, id, rec_chain, pep_chain, rank, model_no]

        model_file = join(output_dir, id + '.fasta', pred_fn)
        native_file = join(input_dir, id + '.pdb')

        cmd.load(native_file, 'native')
        cmd.load(model_file, model_name)

        print(model_name, rec_chain, pep_chain, model_file)
        cmd.select("native_rec", f'native and chain {rec_chain}')
        cmd.select("native_pep", f'native and chain {pep_chain}')

        cmd.select("afold_rec", f'{model_name} and chain {rec_chain}')
        cmd.select("afold_pep", f'{model_name} and chain {pep_chain}')

        # align peptide chains
        super_alignment_pep = cmd.super('afold_pep', 'native_pep')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_pep])

        seq_alignment_pep = cmd.align('afold_pep', 'native_pep')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_pep])

        # align peptide chains backbone
        super_alignment_pep = cmd.super('afold_pep and backbone', 'native_pep and backbone')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_pep])

        seq_alignment_pep = cmd.align('afold_pep and backbone', 'native_pep and backbone')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_pep])

        # super receptor chains
        super_alignment_rec = cmd.super('afold_rec', 'native_rec')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_rec])

        # save the superimposed structure
        cmd.select('model_to_save', model_name)
        super_filename=f'{model_name}_superimposed.pdb'
        print(super_filename)
        save_to_file = os.path.join(input_dir, super_filename)
        cmd.save(save_to_file, model_name, format='pdb')

	# super receptor chain backbones
        super_alignment_rec = cmd.super('afold_rec and backbone', 'native_rec and backbone')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_rec])

        # calculate rmsd-s
        seq_alignment_rec = cmd.align('afold_rec', 'native_rec')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_rec])

	# calculate rmsd by backbone
        seq_alignment_rec = cmd.align('afold_rec and backbone', 'native_rec and backbone')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_rec])

        super_complex = cmd.super(model_name, 'native')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_complex])

        super_complex = cmd.super(f'{model_name} and backbone', 'native and backbone')
        metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_complex])

        #cmd.color('brown', model_name)
        # color receptor interface of model in yellow
        #cmd.color('yellow', 'interface_rec_afold')

	# color peptide of model in cyan
        #cmd.color('cyan', 'afold_pep')
        #cmd.show(representation='sticks', selection='afold_pep')

        metrics.append(metrics_for_a_model)

    #cmd.set_name('native', complex)
    #cmd.save(f'{input_dir}/{complex}.pse', format='pse')

    # create column names
    colnames = ['model_name', 'pdb_id', 'rec_chain', 'pep_chain', 'rank', 'model_no']
    colnames_for_aln = ['rms_after_ref', 'no_aln_atoms_after_ref', 'ref_cycles',
                        'rms_before_ref', 'no_aln_atoms_before_ref', 'raw_score',
                        'no_aln_residues']

    for type in ['_super_pep', '_align_seq_pep', '_super_pep_bb', '_align_seq_pep_bb',
                 '_super_rec', '_super_rec_bb', '_align_seq_rec', '_align_seq_rec_bb',
                 '_complex', '_complex_bb']:

        new_colnames = [s + type for s in colnames_for_aln]
        colnames = colnames + new_colnames

    # saving calculated metrics
    out_csv_dir = os.path.dirname(input_dir.rstrip('/'))
    metrics_df = pandas.DataFrame(metrics, columns = colnames)
    metrics_df.to_csv(os.path.join(out_csv_dir, 'metric' + '.csv'))
