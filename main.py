
import module
import argparse
import configargparse

from parser import parse_args

def prune_extra_atoms(args):
    ''' Remove extra atoms from pdb file.
        Currently only used on native pdb.
    '''
    pipeline = module.pipeline('init_atom_prune', args)
    pipeline.prune_extra_atoms()

def process_input(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''
    pipeline = module.pipeline('init_input_procs', args)
    pipeline.process_input()


def process_input_from_fasta(args):
    ''' Combine multi chains with poly-g linker
          and merge multiple fasta files into one.
        Also store residue id where each single chain starts
          and chain name of each protein
    '''
    pipeline = module.pipeline('init_input_procs_fasta', args)
    pipeline.process_input_from_fasta()


def process_output(args):
    ''' Remove linker from predicted pdb and reorder
          residue number as per gt pdb
        Calculate rmsd and export to csv file
    '''
    pipeline = module.pipeline('init_output_procs', args)
    pipeline.process_output()

def locate_extra_atoms(args):
    ''' Locate extra atoms in gt pdb files '''
    pipeline = module.pipeline('init_atom_locating', args)
    pipeline.locate_extra_atoms()

if __name__ == "__main__":

    parser = configargparse.ArgumentParser()
    config = parse_args(parser)
    args = argparse.Namespace(**config)

    for operation in args.operations:
        if operation == 'process_input':
            if args.from_fasta:
                process_input_from_fasta(args)
            else:
                process_input(args)
        elif operation == 'process_output':
            process_output(args)
        elif operation == 'locate_extra_atoms':
            locate_extra_atoms(args)
        elif operation == 'prune_extra_atoms':
            prune_extra_atoms(args)
        elif operation == 'assert_fasta':
            assert_fasta(args)
