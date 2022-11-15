#!/bin/bash

ds=$1
dir="../../data/bioinfo/output/${ds}_af_full/poly_g_20"
out_s_fn="../../data/bioinfo/input/seq_to_pred/${ds}/tmp_done.txt"

success_pdbs=()

for subdir in $dir/*; do
    if [[ "$subdir" == *"fasta"* ]] ; then
	#echo "${subdir}/features.pkl"
	if [[ -f "${subdir}/features.pkl" ]]; then
	    base_dir=$(basename $subdir)
	    pdb_id=${base_dir:0:4}
	    success_pdbs+=($pdb_id)
	fi
    fi
done

#echo $success_pdbs
printf "%s," "${success_pdbs[@]}" > $out_s_fn
