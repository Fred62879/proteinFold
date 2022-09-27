#!/bin/bash

server=graham
#seqs=(1YCQ.fasta 2AZE.fasta)
#seqs=(1YCQ.fasta 2AZE.fasta 2M3M.fasta 2QTV.fasta 2RSN.fasta 4U7T.fasta)
#seqs=(1AWR.fasta 1CZY.fasta 1EG4.fasta 1ELW.fasta 1ER8.fasta 1JD5.fasta)
#seqs=(1ozs 1tce 1i8h 1h8b 1h3h 1sb0 2laj 2xze 2m3m 2czy 2rsn 3ql9 4x34 4z0r)
#seqs=(1h3h 1h8b 1i8h 1ozs 1sb0 1tce 2czy 2laj 2m3m 2rsn)
#seqs=(3ixs 3fdt 3g2v 5urn 4gng 5heb 3mxc 3wgx 3g2u)
seqs=(3g2u 3mxc 3g2v 3fdt 5urn)

local_data_dir="../../data"
#local_data_dir="/media/fred/Local Disk/Projects/bioinfo/data"
cc_data_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo"

#############
#upload data
#############
upload_data=false

if $upload_data; then
    exps=(idr/polg_g_20)

    for exp in ${exps[@]} ; do
        for pdb_id in ${seqs[@]} ; do
            cc_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo/input/seq_to_pred/${exp}"
            local_dir="../../data/input/${exp}/${pdb_id}"
        done
    done

    scp $local_dir $cc_dir
fi

#################
# download input
#################
download_input=true

cc_input_dir="${cc_data_dir}/input/seq_to_pred"

if $download_input; then
    for exp in ${exps[@]}
    do
        cur_local_data_dir="${local_data_dir}/input/${exp}"

        for fn in pdb_cho0.npy
        do
            cc_fn="${cc_input_dir}/${exp}/${fn}"
            if [ ! -d $local_data_dir ]
            then
                mkdir -p $local_data_dir
            fi

            echo $cc_fn
            echo $cur_local_data_dir
            scp $cc_fn $cur_local_data_dir
        done
    done

    cd "$cur_dir"
fi

##################
# download output
##################
download_files=true

#exps=(idr_af_full/poly_g_96)
exps=(dibs_af_full/poly_g_6 dibs_af_full/poly_g_24 dibs_af_full/poly_g_48 dibs_af_full/poly_g_96)

cc_output_dir="${cc_data_dir}/output"
local_output_dir="${local_data_dir}/output"

if [ ! -d $local_output_dir ]
then
    mkdir -p $local_output_dir
fi

for exp in ${exps[@]}
do
    for seq in ${seqs[@]}
    do
        local_pdb_dir="${local_output_dir}/${exp}/${seq^^}.fasta"
        if [ ! -d $local_pdb_dir ] ; then
            mkdir -p $local_pdb_dir
        fi

        if $download_files ; then
            for nm in ranked_0.pdb ranking_debug.json
            do
                ccdir="${cc_output_dir}/${exp}/${seq^^}.fasta/${nm}"
                #echo $ccdir
                #echo $local_pdb_dir
                scp $ccdir $local_pdb_dir
            done
        fi
    done
done
