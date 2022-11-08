# This is the repo for the Evaluation pipeline of alphafold prediction on idr proteins

## Install
Create a conda virtual environment via `conda env create --file environment.yml`, then activate using `conda activate protein_env`.

## Data Directory
### The pipeline defaults to certain data directory structure shown as follows
<pre>
data_dir
   |
   |_input
   |   |
   |   |_dataset_name
   |   |        |
   |   |        |_pdbs (pdb files can be automatically downloaded)
   |   |        |
   |   |        |_source_fasta (fasta files require manual download)
   |   |
   |   |_dataset_spec (json file)
   |
   |
   |_output
       |
       |_dataset/alphafold_name
                 |
                 |_strategy_name (e.g. poly_g_20)
                        |
                        |_individual pdb prediction (e.g. 1C8O.fasta)
                                     |
                                     |_ranked_0.json (requires manual download)
                                     |
                                     |_ranking_debug.json (requires manual download)
</pre>

### idr dataset 1
This dataset has proteins posted after 2021.
`sampled_data/input/ds1/pdbs` stores three sample pdb proteins. `sampled_data/input/idr/poly_g_20_fasta` stores the corresponding amino acids sequence (fasta format) for each of these three pdb file (poly-g liner length 20, using fasta sequence).

The predictions (ranked_0.pdb) of alphafold are stored in `output/ds1_af_full/poly_g_20_fasta/{pdb_id}.fasta`.

## Run
### Input processing
Input processing deals with pdb and fasta file download and also gathers necessary info about the pdb (chain id, boundary residue id of each chain, and start residue id of each chain for each protein) which are necessary for the later output processing.

#### pdb download and processing
For ground-truth pdbs, either they are pre-downloaded and stored under `data/input/{dataset_name}/pdbs` <sup>1</sup> or automatic download should be enabled. In case of automatic download, if a dataset spec file (json file that specifis, for each pdb id, relevant info incl. chain id etc.) is provided <sup>2</sup>, the pipeline will parse this file to get a list of pdb ids to download. Otherwise <sup>3</sup>, a plain text file of comma separated pdb ids must be provided and stored under `data/input/{dataset_name}` for downloading. Note, sometimes pdb downloading may be interrupted and error will be incurred. If this happened, simply re-run the input processing step. The pipeline will only download the rest pdbs.

<sup>1</sup>: `download_pdbs=False` <br />
<sup>2</sup>: `has_dataset_spec=True, download_pdbs=True` <br />
<sup>3</sup>: `has_dataset_spec=False, download_pdbs=True` <br />

#### fasta download and processing
Fasta download can be skipped if you are using multimer strategy and have already downloaded fasta files for prediction<sup>1</sup>. Otherwise, the pipline supports either generating fasta files via extracting sequences from ground-truth pdb files <sup>2</sup> or users have to manually download fasta files from rcsb.org and store them under `data/input/{dataset_name}/source_fasta`. Note that rcbs.org have separate fasta files for each chain of a single protein. Depending on the prediction strategy, further processing may be necessary. If poly glycine linker strategy is used <sup>3</sup>, two chains of a single protein needs to be linked together by poly glycine residues. If alphafold multimer strategy is used (and you haven't downloaded fasta files yet) <sup>4</sup>, no further processing is required.

<sup>1</sup>: `skip_fasta_download=True, strategy=multimer` <br \>
<sup>2</sup>: `generate_fasta_from_pdb=True` <br \>
<sup>3</sup>: `strategy=poly_g_link` <br \>
<sup>4</sup>: `strategy=multimer` <br \>

### gather other info
For poly glycine linker strategy, as the linker needs to be removed during output processing, we need to record id of the 1st residue for each chain of each protein. This is done during the linking process. If you choose multimer strategy, this processing will not be incurred.

As sequences from downloaded fasta files generally contain a few more residues than the ground-truth pdbs, it's better (required for pymol) that these extra residues are removed<sup>1</sup>. Moreover, alphafold predictions have residue ids being contiguous and starting from 1. While first residue ids of ground-truth pdbs are not necessarily 1 and residue ids between chains may not be contiguous. To renumber alphafold prediction residue ids <sup>2</sup> s.t. they are the same as per ground-truth pdbs, we need to record id of the first and last residue of each chain for each ground-truth pdb.

<sup>1</sup>: `remove_extra_residues=True` <br />
<sup>2</sup>: `renumber_residues=True` <br />


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
