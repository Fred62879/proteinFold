#!/bin/bash

server=beluga
exps=(idr_84_pf_full/poly_g_6)
seqs=(1YCQ.fasta 2AZE.fasta 2M3M.fasta 2QTV.fasta 2RSN.fasta 4U7T.fasta)

#############
#upload data
#############
upload_data=true
cc_dir='fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo/input/seq_to_pred/peptide/poly_g_6'
local_dir='../data/input/peptide/poly_g_6/*'

if $upload_data; then
    scp $local_dir $cc_dir
fi

#################
# download input
#################
download_input=false
input_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo/input"

if $download_input; then
    cur_dir=$(pwd)

    cc_input_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/astro/${img_type}_input"
    input_dir="../../data/${img_type}_input/"
    cd $input_dir

    for dir in "sampled_pixl_ids/cutout_981316_512_0_0_mask_rand_diff_10_279_0"
    do
        ccfn="${cc_input_dir}/${dir}/*"
        if [ ! -d $dir ]
        then
            mkdir -p $dir
        fi
        scp $ccfn $dir
    done

    cd "$cur_dir"
fi

##################
# download output
##################
output_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo/output"

# create output parent dir
data_dir="../data/output/${exp}"

if [ ! -d $data_dir ]
then
    mkdir -p $data_dir
fi
cd $data_dir

download_files=true

for exp in ${exps[@]}
do
    for seq in ${seqs[@]}
    do
        local_out="${exp}/${seq}"
        if [ ! -d $local_out ] ; then
            mkdir -p $local_out
        fi

        if $download_files ; then
            for nm in ranked_0.pdb ranking_debug.json
            do
                ccdir="${output_dir}/${local_out}/${nm}"
                echo $ccdir
                echo $local_out
            scp $ccdir $local_out
            done
        fi
    done
done
