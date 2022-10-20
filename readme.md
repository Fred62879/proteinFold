# This is the repo for the Evaluation pipeline of alphafold prediction on idr proteins

## Directory
### idr dataset 1
This dataset has proteins posted after 2021.
`sampled_data/input/ds1/pdbs` stores three sample pdb proteins. `sampled_data/input/idr/poly_g_20_fasta` stores the corresponding amino acids sequence (fasta format) for each of these three pdb file (poly-g liner length 20, using fasta sequence).

The predictions (ranked_0.pdb) of alphafold are stored in `output/ds1_af_full/poly_g_20_fasta/{pdb_id}.fasta`.

## Install
Create a conda virtual environment via `conda env create --file environment.yml`, then activate using `conda activate protein_env`.

## Run
### Input processing
Currently, the pipeline extract amino acid sequences from the ground-truth pdb file. And also stores all chain identifiers for each protein. To be able to remove linker from predictions, input processing also records range of each original chain when they are being linked into a single chain. These files are stored in `input/idr/poly_g_6`.

### Output processing
Please first download alphafold prediction files yourself, the pipeline currently doesn't support automatic downloading. But I have a few sample alphafold prediction files prepared as described in Directory section.

After downloading the prediction files (ranked_0.pdb) from the remote server, output processing removes linker and re-numbers residue id for each predicted pdb file so that their residue sequence and id are exactly the same as gt.

Then the pipeline can calculate metircs in two ways
1) Using pymol to superimpose the receptor of the ground-truth and predicted protein and then calculate rmsd values for the ligand.
2) Using Biopython to calculate dockQ score (adapted from https://github.com/bjornwallner/DockQ)

Set `dockq` to `True` if you want the 2nd and to `False` for the 1st option.

### Running
`configs/sample.ini` contains all necessary configurations to run the pipeline. For example, to run input processing, simply uncomment `operations=[process_input]` (keep other operations commented).

To run the pipeline, call `python main.py --config configs/sample.ini`.
