
import argparse
import configargparse

from pymol import cmd
from parser import parse_args
from modules import input_procs_module, output_procs_module


def prune_extra_atoms(args):
    ''' Remove extra atoms from pdb file.
        Currently only used on native pdb.
    '''
    pipeline = input_procs_module.pipeline('init_atom_prune', args)
    pipeline.prune_extra_atoms()

def process_input(args):
    ''' Combine multi chains with poly-g linker and merge multiple fasta files into one.
        Also store residue id where each single chain starts and chain name of each protein
    '''
    pipeline = input_procs_module.pipeline('init_input_procs', args)
    pipeline.process_input()

def process_output(args):
    ''' Remove linker from predicted pdb and reorder residue number as per gt pdb
        Calculate rmsd and export to csv file
    '''
    pipeline = output_procs_module.pipeline('init_output_procs', args)
    pipeline.process_output()

def locate_extra_atoms(args):
    ''' Locate extra atoms in gt pdb files '''
    pipeline = input_procs_module.pipeline('init_atom_locating', args)
    pipeline.locate_extra_atoms()

def plot_metrics(args):
    ''' plot metrics against specified variables '''
    pipeline = input_procs_module.pipeline('init_metric_plotting', args)
    pipeline.plot_metrics()


if __name__ == "__main__":
    parser = configargparse.ArgumentParser()
    config = parse_args(parser)
    args = argparse.Namespace(**config)

    for operation in args.operations:
        if operation == 'process_input':
            process_input(args)
        elif operation == 'process_output':
            process_output(args)
        elif operation == 'locate_extra_atoms':
            locate_extra_atoms(args)
        elif operation == 'prune_extra_atoms':
            prune_extra_atoms(args)
        elif operation == 'plot_metrics':
            plot_metrics(args)

    if not args.dockq:
        cmd.quit()
