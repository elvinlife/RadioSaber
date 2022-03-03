#
run_onecase() {
    slice=$1
    for i in $(seq 0 2); do
        let "allue=4*$slice"
        #../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 7 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/nvs_pf$i.log > /dev/null & #${slice}slices/nvs_pf$i.error
        #../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 8 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/greedy_pf$i.log > /dev/null &
        #../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 9 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/subopt_pf$i.log > /dev/null

        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 10 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/maxcell_pf$i.log > ${slice}slices/maxcell_pf$i.error &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 11 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/vogel_pf$i.log > /dev/null
    done
}

run_onecase 4 &
run_onecase 10 &
run_onecase 20 &
#run_onecase 40
#run_onecase 100 &
