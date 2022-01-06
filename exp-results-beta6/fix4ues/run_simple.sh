#
run_onecase() {
    slice=$1
    for i in $(seq 0 2); do
        let "allue=4*$slice"
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 7 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/nvs_pf$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 8 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/greedy_pf$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 9 1 30 0.1 128 $i ${slice}slices/config 2> ${slice}slices/subopt_pf$i.log > /dev/null &
    done
}

#for slice in "4" "10" "20" "40" "100"; do
for slice in "10"; do
  run_onecase $slice
done

