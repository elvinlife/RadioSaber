# fix to 100 flows
run_onecase() {
  s=$1
  for i in $(seq 0 0); do
    ../../LTE-Sim SingleCellCustomize 3 1 0 0 100 0 0 7 1 30 0.1 foreman_H264_128k.dat $i ${s}slices/config 2> ${s}slices/nvs_pf$i.log > /dev/null &
    ../../LTE-Sim SingleCellCustomize 3 1 0 0 100 0 0 8 1 30 0.1 foreman_H264_128k.dat $i ${s}slices/config 2> ${s}slices/greedy_pf$i.log > ${s}slices/greedy_pf$i.error &
  done
}

#run_onecase 50
run_onecase 25
run_onecase 10
