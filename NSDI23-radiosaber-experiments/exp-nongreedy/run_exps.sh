#!/bin/bash
ODIR="ues"

ueArray=("10" "20" "30")

run_onecase() {
    i=$1
    for ue in ${ueArray[@]}; do
        ../../LTE-Sim SingleCellWithI 1 7 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/nvs_${2}$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 9 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/maxcell_${2}$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 11 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/nvs_optimal_${2}$i.log > /dev/null &
    done
}

TRACE_PATH="../../cqi-traces-noise0/"
cp $TRACE_PATH/mapping0.config $TRACE_PATH/mapping.config
run_onecase 2 pf
