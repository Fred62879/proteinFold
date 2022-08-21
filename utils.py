
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(fn):
    fasta_sequences = SeqIO.parse(open(fn),'fasta')
    for fasta in fasta_sequences:
        monomer = str(fasta.seq)
        break
    return monomer

def read_fasta_from_pdb(fn, print=False):
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

def save_fasta(seq, fn, id='0'):
    seq = SeqRecord(Seq(seq),id=id)
    SeqIO.write(seq, fn, 'fasta')

def parse_to_groups(dir):
    prev_id = ''
    cur_group, groups = [], []

    fns = sorted(os.listdir(dir))
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

def poly_g_link(indir, outdir, fasta_groups, n_g):
    poly_g = 'G' * n_g
    chain_start_ids = {}

    for fasta_group in fasta_groups:
        multimer = ''
        cur_chain_start_ids = []
        pdb_id = (fasta_group[0]).split('_')[0]

        for fn in fasta_group:
            monomer = read_fasta(join(indir, fn + '.fasta'))
            cur_chain_start_ids.append(len(multimer))
            multimer += monomer + poly_g

        # add last start_id for convenience of polyg removal (for loop)
        cur_chain_start_ids.append(len(multimer))
        multimer = multimer[:-n_g]

        ofn = join(outdir, pdb_id + '.fasta')
        save_fasta(multimer, ofn, id=pdb_id)

        chain_start_ids[pdb_id] = cur_chain_start_ids
        ofn = join(outdir, 'chain_start_ids.pkl')
        with open(ofn, 'wb') as fp:
            pickle.dump(chain_start_ids, fp)

def remove_linker_rename_chain(pred_fn, out_fn, chain_start_ids, linker_len, chain_names, shift=True):
    ''' need to shift resid id by 1
        otherwise, chain id modification leads to wrong results
    '''
    fasta = read_fasta_from_pdb(pred_fn)[0]
    for i in range(1, len(chain_start_ids)-1):
        id = chain_start_ids[i]
        assert(fasta[id-linker_len : id] == 'GGGGGG')

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pred_fn)
    df = ppdb.df['ATOM']
    linker_atom_ids = []

    acc_num_removed_resids = 0
    n = len(chain_start_ids)

    for i in range(n-1):
        resid_lo = chain_start_ids[i]

        # resid_hi is residue id of first 'g' is fasta
        resid_hi = chain_start_ids[i+1] - linker_len
        if shift: resid_lo += 1;resid_hi += 1
        atom_ids = df[ (df['residue_number'] >= resid_lo) &
                       (df['residue_number'] < resid_hi) ].index
        atom_lo, atom_hi = atom_ids[0], atom_ids[-1]
        print(resid_lo, resid_hi, atom_lo, atom_hi)

        # df.loc upper bound is inclusive
        # (atom_hi is id of last atom of current chain
        df.loc[atom_lo:atom_hi ,'chain_id'] = chain_names[i]
        df.loc[atom_lo:atom_hi ,'residue_number'] -= acc_num_removed_resids

        if i != n - 2:
            acc_num_removed_resids += linker_len

            linker_atom_lo = atom_hi+1
            linker_resid_hi = chain_start_ids[i+1]
            if shift: linker_resid_hi += 1
            linker_atom_hi = df[ (df['residue_number'] < linker_resid_hi) ].index[-1]
            linker_atom_ids += list(np.arange(linker_atom_lo, linker_atom_hi+1))
    '''
    linker_ids = reduce(
        lambda acc, linker_start_id:
        acc + [i for i in range(linker_start_id, linker_start_id+linker_len)],
        linker_start_ids, [])
    '''

    df.drop(linker_atom_ids, axis=0, inplace=True)
    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def split_chains_wo_env(fn, model='native', create=False, select=False, split_only=True, prefix=False):

    cmd.load(fn, model)
    count, chains = 0, []
    models = cmd.get_object_list('(all)')

    for chain in cmd.get_chains('all and model %s' % (model)):
        if chain == '': chain = "''"
        count += 1
        if not prefix: name = f'{model}_{chain}'
        else:          name = '%s%04d' % (prefix, count)

        if create: cmd.create(name, f'{selection} and model {model} and chain {chain}')
        if select: cmd.select(name, f'{selection} and model {model} and chain {chain}')

        chains.append(chain)
        if split_only: cmd.delete(model)
        else:          cmd.disable(model)

    return chains

def split_chains_w_env(selection='(all)', model='native', create=False, select=False, split_only=True, prefix=False):
    count, chains = 0, []
    models = cmd.get_object_list('(' + selection + ')')

    for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
        if chain == '': chain = "''"
        count += 1
        if not prefix: name = f'{model}_{chain}'
        else:          name = '%s%04d' % (prefix, count)

        if create: cmd.create(name, f'{selection} and model {model} and chain {chain}')
        if select: cmd.select(name, f'{selection} and model {model} and chain {chain}')

        chains.append(chain)
        if split_only: cmd.delete(model)
        else:          cmd.disable(model)

    return chains

def find_subseq_range(seq1, seq2):
    n, m = len(seq1), len(seq2)
    for i in range(n):
        if seq1[i:i+m] == seq2:
            return np.array([[0,i],[i+m,n]])
    assert(False)

def find_prune_ranges_all_chains(seqs1, seqs2):
    ''' assume seq in seqs1 is longer than corresponding seq in seq2
    '''
    acc_len = 0
    rnge = []
    for seq1, seq2 in zip(seqs1, seqs2):
        cur_rnge = find_subseq_range(seq1, seq2)
        cur_rnge += acc_len
        rnge.extend(list(cur_rnge))
        acc_len += len(seq1)
    return np.array(rnge)

def prune_seq_given_range(df, rnge, to_remove_atom_ids, shift=True):
    ''' remove atoms falling within residue id lo (inclusive) and hi (exclusive)
    '''
    resid_lo,resid_hi = rnge
    if shift: resid_lo +=1; resid_hi += 1
    atom_ids = df[ (df['residue_number'] >= resid_lo) &
                   (df['residue_number'] < resid_hi) ].index
    to_remove_atom_ids.extend(atom_ids)

def prune_seq_given_ranges(in_fn, out_fn, ranges, shift=True):
    #if len(ranges) == 0: return
    ppdb_pred = PandasPdb()
    _ = ppdb.read_pdb(in_fn)
    df = ppdb.df['ATOM']
    to_remove_atom_ids = []

    for rnge in ranges:
        prune_seq_given_range(df, rnge, to_remove_atom_ids, shift=shift)
    df.drop(to_remove_atom_ids, axis=0, inplace=True)

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def match_seq(fasta1, fasta2, ranges):
    ''' prune fasta2 so it match fasta1 exactly
        fasta1 is a continuous subseq of fasta2
    '''
    resid_lo, resid_hi = find_subseq_range(fasta1, fasta2)

def align_chains(gt_fn, pred_fn, prune_ranges=None):
    ''' prune predicted and gt protein so they have same aa seq
        only used when downloaded fasta doesn't match
        fasta extracted from downloaded pdb
    '''
    if prune_ranges is not None:
        prune_seq_given_ranges(gt_fn, prune_ranges[0])
        prune_seq_given_ranges(pred_fn, prune_ranges[1])
    else:
        assert(False)
        fastas_gt = read_fasta_from_pdb(gt_fn)
        fastas_pred = read_fasta_from_pdb(pred_fn)
        assert(len(fastas_gt) == len(fastas_prede))
        for fasta_gt, fasta_pred in zip(fastas_gt, fastas_pred):
            if len(fasta_gt) > len(fasta_pred):
                match_seq(fasta_gt, fasta_pred)
            elif len(fasta_gt) < len(fasta_pred):
                match_seq(fasta_pred, fasta_gt)
