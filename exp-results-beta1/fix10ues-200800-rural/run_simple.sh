#!/bin/bash
run_onecase() {
    slice=$1
    for i in $(seq 1 2); do
        let "allue=10*$slice"
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 7 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/nvs_pf$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 8 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/greedy_pf$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 9 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/subopt_pf$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 10 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/upperbound_pf$i.log > /dev/null &

        #../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 1 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/single_pf$i.log > /dev/null &
        #../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 8 1 30 0.1 128 $i ${slice}slices/config 2> /dev/null > ${slice}slices/greedy_pf$i_100.error &
    done
}

run_onecase 4 &
run_onecase 10 &
run_onecase 20 &
#run_onecase 32
