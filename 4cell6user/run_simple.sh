for i in $(seq 0 0); do
    timeout 20 ../LTE-Sim SingleCellWithI 4 1 6 0 0 1 0 7 1 30 0.1 128 $i 2> nvs_tta$i.log > /dev/null
    timeout 30 ../LTE-Sim SingleCellWithI 4 1 6 0 0 1 0 8 1 30 0.1 128 $i 2> oracle_tta$i.log > /dev/null
done
