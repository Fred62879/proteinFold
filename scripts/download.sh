#!/bin/bash

server=graham
ds=ds2
seqs=(5NCU 2QME 2FOM 4IF6 1GL0 4UX6 4FDD 2HKQ 2APO 2RVQ 2F91 3O40 2XPP 3O43 4GF3 3B0Z 2F31 2BZW 2VDA 5GNA 3MP7 3PLV 2C5I 5WUJ 1WDG 5WRV 1B0N 3UE5 2L27 1P4Q 2NNU 2MLZ 2LE8 1MIU 5NQF 3S1B 1MZW 1M46 3ZGA 2MNJ 3FDL 5IKJ 2BB3 2L9S 2VOH 5F5S 2XJY 4DZU 2KXW 1SYQ 2N4Q 2KXQ 1WKW 5CHL 2RIV 1O9A 5OWU 4JHK 2K9U 1SGH 3O42 1SP4 3PLX 2RAW 3C0T 2RR3 1MTP 1Y43 3GME 2MLX 2YLE 5IXD 1A2X 3EAR 2N01 5GK9 5U4K 2LFW 2LYD 1UHD 2K2I 1YYP 1ETR 2M8S 5ABX 2JXC 1HLE 5VMO 1JDH 6F9S 2Y0M 3DBO 2L2I 3E5A 4DZV 5F74 5YI8 2LPB 2KT5 1MQS 2A7U 2VU8 3ZLD 4IKA 2RQW 5XIU 2L1C 5H65 9PAI 1G0Y 5JW9 4D2L 1H8B 3N00 2KWK 3FII 1MCT 3QRX 2V6X 2DYO 3ZWZ 5YAH 4D0G 2L53 6F9W 3AJB 1MK2 5FCG 2XBX 4ZRL 2MFQ 1QGK 3AU4 1RF3 5CUL 3TZ1 3GP2 2LXM 2L8T 1BH8 2FID 5H7Y 4BRU 3KL4 1VC3 2LUH 3BTP 1MCV 2L9B 2RDD 2IPP 4CXF 2BTC 4DXR 2L1W 5WN9 1HTR 1SC5 2HQW 2YBF 3NDF 5KIY 2LI5 2CO4 4C5H 4CLQ 4NUT 1QTX 1S6C 5GTB 2GDC)

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
download_input=false

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
            scp $cc_fn $cur_local_data_dir
        done
    done

    cd "$cur_dir"
fi

##################
# download output
##################
download_files=true

exps=(${ds}_af_full/poly_g_20_fasta)

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
        echo $seq
        local_pdb_dir="${local_output_dir}/${exp}/${seq^^}.fasta"
        if [ ! -d $local_pdb_dir ] ; then
            mkdir -p $local_pdb_dir
        fi

        if $download_files ; then
            for nm in ranked_0.pdb ranking_debug.json
            do
                ccdir="${cc_output_dir}/${exp}/${seq^^}.fasta/${nm}"
                scp $ccdir $local_pdb_dir
            done
        fi
    done
done
