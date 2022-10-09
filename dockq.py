
import os
import re
import Bio.PDB
import numpy as np

from Bio.SVDSuperimposer import SVDSuperimposer


def parse_fnat(fnat_out):
    inter = []
    fnat, fnonnat, nonnat_count = -1, -1, -1
    nat_total, model_total, nat_correct = -1, -1, -1

    for line in fnat_out.split("\n"):
        line = line.rstrip('\n')
        match = re.search(r'NATIVE: (\d+)(\w) (\d+)(\w)',line)

        if re.search(r'^Fnat',line):
            list = line.split(' ')
            fnat = float(list[3])
            nat_correct = int(list[1])
            nat_total = int(list[2])
        elif re.search(r'^Fnonnat',line):
            list = line.split(' ')
            fnonnat = float(list[3])
            nonnat_count = int(list[1])
            model_total = int(list[2])
        elif match:
            res1, chain1 = match.group(1), match.group(2)
            res2, chain2 = match.group(3), match.group(4)
            inter.append(res1 + chain1)
            inter.append(res2 + chain2)

    return (fnat, nat_correct, nat_total, fnonnat, nonnat_count, model_total, inter)

def get_fnat(model, native, exec_path, capri_peptide=False):
    cmd_fnat = exec_path + '/fnat ' + model + ' ' + native + ' 5 -all'
    cmd_interface = exec_path + '/fnat ' + model + ' ' + native + ' 10 -all'

    if capri_peptide:
        cmd_fnat = exec_path + '/fnat ' + model + ' ' + native + ' 4 -all'
        cmd_interface = exec_path + '/fnat ' + model + ' ' + native + ' 8 -cb'

    fnat_out = os.popen(cmd_fnat).read()

    (fnat, nat_correct, nat_total, fnonnat, nonnat_count, model_total, interface5A) = parse_fnat(fnat_out)
    assert fnat != -1, "Error running cmd: %s\n" % (cmd_fnat)

    inter_out = os.popen(cmd_interface).read()
    (fnat_bb,nat_correct_bb,nat_total_bb,fnonnat_bb,nonnat_count_bb,model_total_bb,interface) = parse_fnat(inter_out)
    assert fnat_bb != -1, "Error running cmd: %s\n" % (cmd_interface)

    info = { 'fnat': fnat,
             'fnonnat': fnonnat,
             'nat_total': nat_total,
             'nat_correct': nat_correct,
             'model_total': model_total,
             'nonnat_count': nonnat_count,
             'interface': interface
             }
    #return interface5A,fnat_bb,nat_correct_bb,nat_total_bb,fnonnat_bb,nonnat_count_bb,model_total_bb
    return info

def get_sample_atoms(atom_for_sup, sample_model):
    atoms_def_sample=[]
    for sample_chain in sample_model:
        chain = sample_chain.id
        for sample_res in sample_chain:
            is_het_atm = sample_res.get_id()[0] != ' '
            if is_het_atm: continue

            resname = sample_res.get_id()[1]
            key = str(resname) + chain
            for a in atom_for_sup:
                atom_key = key + '.' + a
                if a in sample_res:
                    atoms_def_sample.append(atom_key)
    return atoms_def_sample

def get_ref_sample_common_atoms(atom_for_sup, atoms_def_sample, ref_model):
    atoms_def_in_both = []
    for ref_chain in ref_model:
        chain = ref_chain.id
        for ref_res in ref_chain:
            is_het_atm = ref_res.get_id()[0] != ' '
            if is_het_atm:
                # print ref_res.get_id()
                continue

            resname=ref_res.get_id()[1]
            key=str(resname) + chain
            for a in atom_for_sup:
                atom_key=key + '.' + a
                if a in ref_res and atom_key in atoms_def_sample:
                    if atom_key in atoms_def_in_both:
                        print(atom_key + ' already added (Native)!!!')
                    atoms_def_in_both.append(atom_key)
    return atoms_def_in_both

def get_sample_interface(atom_for_sup, sample_model, interface, atoms_def_in_both):
    ''' Return: chain_res - residue for each chain of sample pdb
                sample_atoms - atoms for interface residues in sample pdb
                common_interface - interface residues
    '''
    chain_res, sample_atoms, common_interface = {}, [], []

    for sample_chain in sample_model:
        chain = sample_chain.id
        if chain not in list(chain_res.keys()):
            chain_res[chain]=[]

        for sample_res in sample_chain:
            is_het_atm = sample_res.get_id()[0] != ' '
            if is_het_atm: continue

            resname = sample_res.get_id()[1]
            key = str(resname) + chain
            chain_res[chain].append(key)

            if key in interface:
                for a in atom_for_sup:
                    atom_key = key + '.' + a
                    if a in sample_res and atom_key in atoms_def_in_both:
                        sample_atoms.append(sample_res[a])
                common_interface.append(key)

    return chain_res, sample_atoms, common_interface

def get_ref_interface(ref_model, sample_res, atom_for_sup, atoms_def_in_both, common_interface):
    ''' Return: chain_ref -
                ref_atoms - atoms for interface residues in ref pdb
                common_residues - residues in both sample and ref
    '''
    chain_ref, ref_atoms, common_residues = {}, [], []

    for ref_chain in ref_model:
        chain = ref_chain.id
        if chain not in list(chain_ref.keys()):
            chain_ref[chain] = []

        for ref_res in ref_chain:
            is_het_atm = ref_res.get_id()[0] != ' '
            if is_het_atm: continue

            resname = ref_res.get_id()[1]
            res_key = str(resname) + chain

            res_in_sample = res_key in sample_res[chain]
            if res_in_sample:
                for a in atom_for_sup:
                    atom_key = res_key + '.' + a
                    if a in ref_res and atom_key in atoms_def_in_both:
                        chain_ref[chain].append(ref_res[a])
                        # ????? indent wrong ????
                common_residues.append(res_key)

            if res_key in common_interface:
                # Check if residue number ( .get_id() ) is in the list
                # Append CA atom to list
                for a in atom_for_sup:
                    atom_key = res_key + '.' + a
                    if a in ref_res and atom_key in atoms_def_in_both:
                        ref_atoms.append(ref_res[a])

    return chain_ref, ref_atoms, set(common_residues)

def get_sample_common(sample_model, atom_for_sup, atoms_def_in_both, common_residues):
    ''' Return: chain_sample - atoms of common residues in sample
    '''
    chain_sample = {}
    for sample_chain in sample_model:
        chain = sample_chain.id
        if chain not in list(chain_sample.keys()):
            chain_sample[chain] = []

        for sample_res in sample_chain:
            is_het_atm = sample_res.get_id()[0] != ' '
            if is_het_atm: continue

            resname = sample_res.get_id()[1]
            key = str(resname) + chain
            if key in common_residues:
                for a in atom_for_sup:
                    atom_key = key + '.' + a
                    if a in sample_res and atom_key in atoms_def_in_both:
                        chain_sample[chain].append(sample_res[a])
    return chain_sample

def assert_selections(receptor_chain, ligand_chain, chain_ref, chain_sample):
    assert len(chain_ref[receptor_chain]) == len(chain_sample[receptor_chain]), \
        "Different number of atoms in native and model receptor (chain %c) %d %d\n" % \
        (receptor_chain, len(chain_ref[receptor_chain]), len(chain_sample[receptor_chain]))

    assert len(chain_ref[ligand_chain]) != 0 or len(chain_sample[ligand_chain]) != 0, \
        "Zero number of equivalent atoms in native and model ligand (chain %s) %d %d.\anCheck that the residue numbers in model and native is consistent\n" % \
        (ligand_chain, len(chain_ref[ligand_chain]), len(chain_sample[ligand_chain]))

    assert len(chain_ref[ligand_chain]) == len(chain_sample[ligand_chain]), \
        "Different number of atoms in native and model ligand (chain %c) %d %d\n" % \
        (ligand_chain, len(chain_ref[ligand_chain]), len(chain_sample[ligand_chain]))

def load_and_select(atom_for_sup, sample_model, ref_model, info):
    interface = info['interface']

    atoms_def_sample = get_sample_atoms(atom_for_sup, sample_model)
    atoms_def_in_both = get_ref_sample_common_atoms\
        (atom_for_sup, atoms_def_sample, ref_model)
    sample_res, sample_atoms, common_interface = get_sample_interface\
        (atom_for_sup, sample_model, interface, atoms_def_in_both)
    chain_ref, ref_atoms, common_residues = get_ref_interface\
        (ref_model, sample_res, atom_for_sup, atoms_def_in_both, common_interface)
    chain_sample = get_sample_common\
        (sample_model, atom_for_sup, atoms_def_in_both, common_residues)

    assert len(ref_atoms) != 0, "length of native is zero"
    assert len(sample_atoms) != 0, "length of model is zero"
    assert len(ref_atoms) == len(sample_atoms), "Different number of atoms in native and model %d %d\n" % (len(ref_atoms), len(sample_atoms))

    # reorder chains to make sure ligand is shorter than receptor
    (chain1, chain2) = list(chain_sample.keys())

    len1 = len(sample_res[chain1])
    len2 = len(sample_res[chain2])
    assert len1 != 0, "%s chain has zero length!\n" % chain1
    assert len2 != 0, "%s chain has zero length!\n" % chain2

    class1, class2 = 'ligand', 'receptor'
    ligand_chain, receptor_chain = chain1, chain2
    if (len(chain_sample[chain1]) > len(chain_sample[chain2])):
        receptor_chain, ligand_chain = chain1, chain2
        class1, class2 = 'receptor', 'ligand'

    assert_selections(receptor_chain, ligand_chain, chain_ref, chain_sample)

    info['len1'] = len1
    info['len2'] = len2
    info['chain1'] = chain1
    info['chain2'] = chain2
    info['class1'] = class1
    info['class2'] = class2
    info['ligand_chain'] = ligand_chain
    info['receptor_chain'] = receptor_chain
    return info, ref_atoms, sample_atoms, chain_ref, chain_sample


def calc_DockQ(info, ref_model, sample_model, ref_atoms, sample_atoms, chain_ref, chain_sample):
    fnat = info['fnat']
    ligand_chain = info['ligand_chain']
    receptor_chain = info['receptor_chain']

    # calculate irms (rmsd between sample and ref for whole protein)
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    irms = super_imposer.rms
    info['irms'] = irms

    # calculate Lrms (align receptor and calculate rmsd for ligand)
    super_imposer.set_atoms(chain_ref[receptor_chain], chain_sample[receptor_chain])
    super_imposer.apply(sample_model.get_atoms())
    receptor_chain_rms = super_imposer.rms

    coord1 = np.array([atom.coord for atom in chain_ref[ligand_chain]])
    coord2 = np.array([atom.coord for atom in chain_sample[ligand_chain]])
    #coord1=np.array([atom.coord for atom in chain_ref[receptor_chain]])
    #coord2=np.array([atom.coord for atom in chain_sample[receptor_chain]])

    sup = SVDSuperimposer()
    Lrms = sup._rms(coord1, coord2) # no superimpose
    info['Lrms'] = Lrms

    #super_imposer.set_atoms(chain_ref[ligand_chain], chain_sample[ligand_chain])
    #super_imposer.apply(sample_model.get_atoms())
    #coord1=np.array([atom.coord for atom in chain_ref[receptor_chain]])
    #coord2=np.array([atom.coord for atom in chain_sample[receptor_chain]])
    #Rrms= sup._rms(coord1,coord2)
    #should give same result as above line
    #diff = coord1-coord2
    #l = len(diff) #number of atoms
    #from math import sqrt
    #print sqrt(sum(sum(diff*diff))/l)
    #print np.sqrt(np.sum(diff**2)/l)

    # calculate dockq
    DockQ=(float(fnat) + 1/(1+(irms/1.5)*(irms/1.5)) + 1/(1+(Lrms/8.5)*(Lrms/8.5)))/3
    info['dockQ'] = DockQ
    return info

def calc_metrics(model, native, exec_path, use_CA_only=False, capri_peptide=False):
    atom_for_sup = ['CA','C','N','O']
    if (use_CA_only): atom_for_sup = ['CA']

    # load structure
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure("reference", native)
    sample_structure = pdb_parser.get_structure("model", model)

    # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]

    # calcaulte fnat
    info = get_fnat(model, native, exec_path)

    # load and select atoms and residues
    info, ref_atoms, sample_atoms, chain_ref, chain_sample = load_and_select\
        (atom_for_sup, sample_model, ref_model, info)

    # calculate irms, Lrms, and dockq
    info = calc_DockQ(info, ref_model, sample_model, ref_atoms, sample_atoms, chain_ref, chain_sample)

    return info
