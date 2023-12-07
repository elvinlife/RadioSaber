#!/bin/bash
# ODIR="ues"
ODIR="ues-ip"

ueArray=("5" "10" "15" "20")

run_onecase() {
    i=$1
    for ue in ${ueArray[@]}; do
        ../../LTE-Sim SingleCellWithI 1 7 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/nvs_${2}$i.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 9 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/maxcell_${2}$i.log > /dev/null &
        # ../../LTE-Sim SingleCellWithI 1 8 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/sequential_${2}$i.log > /dev/null &
        # ../../LTE-Sim SingleCellWithI 1 10 1 30 $i 12 ${ue}${ODIR}/config-${2}.json 2> ${ue}${ODIR}/upperbound_${2}$i.log > /dev/null &
    done
}

TRACE_PATH="../../cqi-traces-noise0/"
cp $TRACE_PATH/mapping0.config $TRACE_PATH/mapping.config
run_onecase 0 pf
