#!/bin/bash

upload_data=false
server_to_server=true
download_input=false
download_output=false

ds=$1
server=graham
to_server=narval
#seqs=(1MIU 1MQS 1QGK 2RDD 2VDA 3AU4 3BTP 3GME 4C5H 4FDD 4H1W 4IKA 4U4C 5IKJ 5IXD 5OWU)
#seqs=(1K9O 3HM6 4PJU 4TSH 4XDN 5BWD 5MZ6 5U1T 5V8K 5YVW 5YVY 6BWB 6KBM 6LRY 6OF7 6RM9 6SH6 6XMU 7E2H 7P9U 7PP2 7SR9 7VSI)
seqs=(3HM6 4FZ1 6KBM 6OF7 6XMU 7DCP 7E2H 7MVY 7P9U 7PP2 7TYR 7VSI)


#############
#upload data
#############

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


###################
# server to server
###################

if $server_to_server; then
    local_dir="/home/fred862/scratch/fred862/data/bioinfo/output"
    out_server_dir="fred862@${to_server}.computecanada.ca:~/scratch/fred862/data/bioinfo/output"

    exps=("${ds}_af_full/poly_g_20")
    for exp in ${exps[@]}
    do
	cd "${local_dir}/${exp}"
	to_zip="to_zip"
	mkdir $to_zip

	for seq in ${seqs[@]}
	do
	    echo $seq
	    local_pdb_dir="${seq^^}.fasta"
	    find $local_pdb_dir -name 'features.pkl' -exec cp --parents \{\} /target \;
	    cp --parent "${local_pdb_dir}/features.pkl" $to_zip/
	done

	zip -r "${to_zip}.zip" "${to_zip}/"
	out_pdb_dir="${out_server_dir}/${exp}"
	scp "${to_zip}.zip" $out_pdb_dir
    done
fi

#################
# download input
#################

if $download_input; then
    cc_data_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo"
    local_data_dir="../../data"
    #local_data_dir="/media/fred/Local Disk/Projects/bioinfo/data"

    cc_input_dir="${cc_data_dir}/input/seq_to_pred"

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
            scp $cc_fn $cur_local_data_dir
        done
    done

    cd "$cur_dir"
fi

##################
# download output
##################

if $download_output ; then
    cc_data_dir="fred862@${server}.computecanada.ca:~/scratch/fred862/data/bioinfo"
    local_data_dir="../../data"
    #local_data_dir="/media/fred/Local Disk/Projects/bioinfo/data"

    for exp in ${exps[@]}
    do
	exps=(${ds}_af_full/poly_g_20_fasta)

	cc_output_dir="${cc_data_dir}/output"
	local_output_dir="${local_data_dir}/output"

	if [ ! -d $local_output_dir ]
	then
	    mkdir -p $local_output_dir
	fi

	for seq in ${seqs[@]}
	do
	    echo $seq
	    local_pdb_dir="${local_output_dir}/${exp}/${seq^^}.fasta"
	    if [ ! -d $local_pdb_dir ] ; then
		mkdir -p $local_pdb_dir
	    fi

	    for nm in ranked_0.pdb ranking_debug.json
	    do
		ccdir="${cc_output_dir}/${exp}/${seq^^}.fasta/${nm}"
		scp $ccdir $local_pdb_dir
	    done
	done
    done

fi
