# This is the repo for the Evaluation pipeline of alphafold prediction on idr proteins

## Directory
`sampled_data/input/idr/pdbs` stores three sample pdb proteins. `sampled_data/input/idr/poly_g_6` stores the corresponding amino acids sequence (fasta format) for each of these three pdb file.

The predictions (ranked_0.pdb) of alphafold are stored in `output/idr_af_full/poly_g_6/{pdb_id}.fasta`.


## Install
Create a conda virtual environment via `conda env create --file environment.yml`, then activate using `conda activate protein_env`.


## Run
### Input processing
Currently, the pipeline extract amino acid sequences from the ground-truth pdb file. And also stores all chain identifiers for each protein. To be able to remove linker from predictions, input processing also records range of each original chain when they are being linked into a single chain. These files are stored in `input/idr/poly_g_6`.

### Output processing
After downloading alphafold predictions (ranked_0.pdb) from the remote server (sample files prepared as described in Directory section), output processing removes linker and re-number residue number from each predicted pdb file. Then the pipeline superimpose the receptor of the ground-truth and predicted protein and then calculate a few rmsd values.

### Running
`configs/test.ini` contains all necessary configurations to run the pipeline. For example, to run input processing, simply do `operations=[process_input]`, and `operations=[process_output]` for output processing.

To run input processing, call `python main.py --config configs/sample.ini`.
