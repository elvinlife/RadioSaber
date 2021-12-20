for ue in "8"; do
#for ue in "2" "4" "6" "8"; do
    for i in $(seq 4 10); do
        let "allue=4*$ue"
        ../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 1 1 30 0.1 128 $i ${ue}ues/config 2> ${ue}ues/single_pf$i.log > /dev/null &
        ../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 7 1 30 0.1 128 $i ${ue}ues/config 2> ${ue}ues/nvs_pf$i.log > /dev/null &
        ../LTE-Sim SingleCellWithI 3 1 $allue 0 0 1 0 8 1 30 0.1 128 $i ${ue}ues/config 2> ${ue}ues/oracle_pf$i.log > /dev/null
    done
done
