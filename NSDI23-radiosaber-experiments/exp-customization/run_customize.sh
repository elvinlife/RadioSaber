#!/bin/bash
ODIR="exp-customize-20slices"

run_onecase() {
    for i in $(seq 0 2); do
        TRACE_PATH="../../cqi-traces-noise0/"
        cp $TRACE_PATH/mapping$i.config $TRACE_PATH/mapping.config

        ../../LTE-Sim SingleCellWithI 1 1 1 30 $i 40 ${ODIR}/config.json 2> ${ODIR}/single_${i}.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 7 1 30 $i 40 ${ODIR}/config.json 2> ${ODIR}/nvs_${i}.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 9 1 30 $i 40 ${ODIR}/config.json 2> ${ODIR}/maxcell_${i}.log > /dev/null &
        sleep 1
    done
}

run_onecase
