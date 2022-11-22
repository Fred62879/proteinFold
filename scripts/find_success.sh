#!/bin/bash

ds=$1
mode=$2
dir="../../../../data/bioinfo/output/${ds}_af_full/poly_g_20"
out_s_fn="../../../../data/bioinfo/input/seq_to_pred/${ds}/tmp_${mode}_done.txt"

success_pdbs=()

for subdir in $dir/*; do
    if [[ "$subdir" == *"fasta"* ]] ; then
	if [[ $mode == "cpu" ]] ; then
	    #echo "${subdir}/features.pkl"
	    if [[ -f "${subdir}/features.pkl" ]]; then
		base_dir=$(basename $subdir)
		pdb_id=${base_dir:0:4}
		success_pdbs+=($pdb_id)
	    fi
	else
	    if [[ -f "${subdir}/ranking_debug.json" && -f "${subdir}/ranked_0.pdb" ]]; then
		echo $subdir
		base_dir=$(basename $subdir)
		pdb_id=${base_dir:0:4}
		success_pdbs+=($pdb_id)
	    fi
	fi
    fi
done

#echo $success_pdbs
printf "%s," "${success_pdbs[@]}" > $out_s_fn
