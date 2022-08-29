
def read_residue_id_from_pdb(pdb_fn):
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(pdb_fn)
    df = ppdb.df['ATOM']
    n = df.shape[0]

    prev_id = df.loc[0,'residue_number']
    resid_ids = [prev_id]

    for atom_id in range(1, n):
       cur_id = df.loc[atom_id, 'residue_number']
       if cur_id != prev_id:
           resid_ids.append(cur_id)
           prev_id = cur_id

    return resid_ids

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

def find_subseq_range(seq1, seq2):
    n, m = len(seq1), len(seq2)
    for i in range(n):
        if seq1[i:i+m] == seq2:
            return np.array([[0,i],[i+m,n]])
    assert(False)

def find_x_range(seq):
    ''' find range of unknown aa subseq in seq '''
    start, in_X, n = -1, False, len(seq)
    ranges = []
    for i in range(n):
        if seq[i] == 'X':
            if not in_X:
                in_X, start = True, i
        elif in_X:
            ranges.append([start,i])
            start, in_X = -1, False

    if start != -1:
        ranges.append([start,n])
    return np.array(ranges)

def find_prune_ranges_all_chains(seqs1, seqs2, prune_X=True):
    ''' assume seq in seqs1 is longer than corresponding seq in seq2
    '''
    acc_len = 0
    rnge1, rnge2 = [], []
    for seq1, seq2 in zip(seqs1, seqs2):
        if prune_X:
            cur_rnge1 = find_x_range(seq1)
            if len(cur_rnge1) != 0:
                rnge1.extend(list(cur_rnge1 + acc_len))
            cur_rnge2 = find_x_range(seq2)
            if len(cur_rnge2) != 0:
                rnge2.extend(list(cur_rnge2 + acc_len))
        else:
            cur_rnge = find_subseq_range(seq1, seq2)
            if len(cur_rnge) != 0:
                rnge1.extend(list(cur_rnge + acc_len))
        acc_len += len(seq1)

    return (np.array(rnge1), np.array(rnge2))

def prune_seq_given_range(df, rnge, to_remove_atom_ids, shift=True):
    ''' remove atoms falling within residue id lo (inclusive) and hi (exclusive)
    '''
    (resid_lo,resid_hi) = rnge
    if resid_lo == resid_hi: return

    if shift: resid_lo +=1; resid_hi += 1
    atom_ids = df[ (df['residue_number'] >= resid_lo) &
                   (df['residue_number'] < resid_hi) ].index
    to_remove_atom_ids.extend(atom_ids)

def prune_seq_given_ranges(in_fn, out_fn, ranges, shift=True):
    #if len(ranges) == 0: return
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fn)
    df = ppdb.df['ATOM']
    to_remove_atom_ids = []

    for rnge in ranges:
        prune_seq_given_range(df, rnge, to_remove_atom_ids, shift=shift)
    df.drop(to_remove_atom_ids, axis=0, inplace=True)

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)
    print(read_fasta_from_pdb(out_fn))

def match_seq(fasta1, fasta2, ranges):
    ''' prune fasta2 so it match fasta1 exactly
        fasta1 is a continuous subseq of fasta2
    '''
    resid_lo, resid_hi = find_subseq_range(fasta1, fasta2)

def align_chains(gt_fn, gt_out_fn, pred_fn, pred_out_fn, prune_ranges=None):
    ''' prune predicted and gt protein so they have same aa seq
        only used when downloaded fasta doesn't match
        fasta extracted from downloaded pdb
    '''
    if prune_ranges is not None:
        print(gt_fn, gt_out_fn)
        print(prune_ranges[1])
        prune_seq_given_ranges(gt_fn, gt_out_fn, prune_ranges[1])
        prune_seq_given_ranges(pred_fn, pred_out_fn, prune_ranges[0])
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

def remove_linker_rename_chain(pred_fn, out_fn, chain_start_ids, linker_len, chain_names):
    ''' rename residue continuously (1-based) between chains
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
        resid_lo += 1;resid_hi += 1
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
            linker_resid_hi += 1
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

def reorder_residue(in_fn, out_fn):
    ''' only used to reorder gt pdb residue id from 1 to n '''

    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fn)
    df = ppdb.df['ATOM']
    n = df.shape[0]

    start, old_id, new_id = 0, df.loc[0,'residue_number'], 1
    for i in range(1, n):
        id = df.loc[i,'residue_number']
        if id != old_id:
            df.loc[start:i,'residue_number'] = new_id
            new_id += 1
            start = i
            old_id = id
    df.loc[start:,'residue_number'] = new_id

    ppdb.df['ATOM'] = df
    ppdb.to_pdb(out_fn)

def reorder_residues(dir, pdb_ids, reorder_fn_str):
    for pdb_id in pdb_ids:
        gt_fn = join(dir, pdb_id + '.pdb')
        out_fn = join(dir, pdb_id + '_' + reorder_fn_str + '.pdb')
        reorder_residue(gt_fn, out_fn)

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

def poly_g_link(indir, outdir, chain_start_ids, fasta_group, poly_g, n_g):
    multimer = ''
    cur_chain_start_ids = []
    pdb_id = (fasta_group[0]).split('_')[0]

    for fn in fasta_group:
        monomer = read_fasta(join(indir, fn + '.fasta'))
        cur_chain_start_ids.append(len(multimer))
        multimer += monomer + poly_g

    # add last start_id for convenience of polyg removal (for loop)
    cur_chain_start_ids.append(len(multimer) + 1) # resid id 1-based
    multimer = multimer[:-n_g]

    ofn = join(outdir, pdb_id + '.fasta')
    save_fasta(multimer, ofn, id=pdb_id)
    chain_start_ids[pdb_id] = cur_chain_start_ids

def poly_g_link_all(indir, outdir, fasta_groups, n_g):
    poly_g = 'G' * n_g
    chain_start_ids = {}

    for fasta_group in fasta_groups:
        poly_g_link(indir, outdir, chain_start_ids, fasta_group, poly_g, n_g)

    return chain_start_ids


#######
# main

def process_output(args):

    with open(args.chain_names_fn, 'rb') as fp:
        chain_names = pickle.load(fp)

    with popen(args.chain_start_ids_fn, 'rb') as fp:
        chain_start_ids = pickle.load(fp)

    with popen(args.gt_chain_bd_resid_ids_fn, 'rb') as fp:
        gt_chain_bd_ids = pickle.load(fp)

    for pdb_id in args.pdb_ids:
        print()
        print(pdb_id)
        dir = join(args.output_dir, pdb_id + '.fasta')
        if not exists(dir):
            print(f'directory {dir} doesn\'t exist')
            continue

        pred_fn = join(dir, args.pred_fn_str)
        gt_pdb_fn = join(args.input_pdb_dir, pdb_id+ '.pdb')
        pred_removed_linker_fn = join(dir, args.removed_linker_fn_str)

        utils.remove_linker(pred_fn, pred_removed_linker_fn, args.n_g,
                            chain_names[pdb_id], chain_start_ids[pdb_id],
                            gt_chain_bd_ids[pdb_id])

        # prune unknown aa from sequences
        align_option = 'pruned_X'
        seqs_pred = utils.read_fasta_from_pdb(pred_removed_linker_fn)
        seqs_gt = utils.read_fasta_from_pdb(gt_pdb_fn)
        print(seqs_pred)
        print(seqs_gt)
        prune_ranges = utils.find_prune_ranges_all_chains(seqs_pred, seqs_gt, prune_X=True)

        #print(prune_ranges)
        gt_removed_X_fn = gt_pdb_fn[:-4] + '_' + align_option + '.pdb'
        pred_removed_X_fn = pred_fn[:-4] + '_' + align_option + '.pdb'
        utils.align_chains(gt_pdb_fn, gt_removed_X_fn, pred_removed_linker_fn, pred_removed_X_fn, prune_ranges)

        # prune extra aa
        align_option = 'aligned'
        seqs_gt = utils.read_fasta_from_pdb(gt_removed_X_fn)
        seqs_pred = utils.read_fasta_from_pdb(pred_removed_X_fn)
        print(seqs_gt)
        print(seqs_pred)

        gt_aligned_fn = gt_pdb_fn[:-4] + '_' + align_option + '.pdb'
        pred_aligned_fn = pred_fn[:-4] + '_' + align_option + '.pdb'
        prune_ranges = utils.find_prune_ranges_all_chains(seqs_pred, seqs_gt, prune_X=False)
        print(prune_ranges)
        utils.align_chains(gt_removed_X_fn, gt_aligned_fn, pred_removed_X_fn, pred_aligned_fn, align_option, prune_ranges)

    cmd.quit()

def process_input(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''

    groups = utils.parse_to_groups(args.source_fasta_dir)
    generate_fasta_from_pdb(in_dir, out_dir, pdb_ids, n_g):

    gt_chain_bd_resid_ids = utils.read_bd_resid_id_all \
        (args.input_pdb_dir, args.pdb_ids)

    with open(args.gt_chain_bd_resid_ids_fn, 'wb') as fp:
        pickle.dump(gt_chain_bd_resid_ids, fp)

    # link multi chains with poly-g link
    chain_start_ids = utils.poly_g_link_all \
        (args.source_fasta_dir, args.linker_fasta_dir, groups, args.n_g)

    # save residue id where each chain start and ends
    with open(args.chain_start_ids_fn, 'wb') as fp:
        pickle.dump(chain_start_ids, fp)

    # get chain names and store locally
    chain_names = utils.get_chain_identifiers_all \
        (args.pdb_ids, args.input_pdb_dir)
    with open(args.chain_names_fn, 'wb') as fp:
        pickle.dump(chain_names, fp)

    #cmd.quit()
