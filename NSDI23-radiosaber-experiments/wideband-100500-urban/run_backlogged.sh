#!/bin/bash
ODIR="randues"

run_onecase() {
    for i in $(seq 0 2); do

        ../../LTE-Sim SingleCellWithI 1 7 1 30 $i 12 ${ODIR}/config.json 2> ${ODIR}/nvs_${i}.log > /dev/null &
        ../../LTE-Sim SingleCellWithI 1 9 1 30 $i 12 ${ODIR}/config.json 2> ${ODIR}/maxcell_${i}.log > /dev/null &
    done
}

run_onecase
