
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

def trim_x(seq):
    ''' remove 'X' from sequence '''
    return seq.replace('X','')

def remove_x_from_pdb(dir, pdb_ids):
    for pdb_id in pdb_ids:
        in_fn = join(dir, pdb_id + '.pdb')
        out_fn = join(dir, pdb_id + '.pdb')
        seq = read_residue_from_pdb(in_fn)
        ranges = find_x_range(seq)
        prune_seq_given_ranges(in_fn, out_fn, ranges)

def find_x_range(seq):
    ''' find range of unknown aa subseq in seq
        output range is lower incluside, upper exclusive
    '''
    start, in_X, n = -1, False, len(seq)
    ranges = []
    for i in range(n):
        if seq[i] == 'X':
            if not in_X:
                in_X, start = True, i+1
        elif in_X:
            ranges.append([start,i+1])
            start, in_X = -1, False

    if start != -1:
        ranges.append([start,n])
    return np.array(ranges)

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
