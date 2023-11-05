#!/bin/bash

TRACE_PATH="../../cqi-traces-noise0"

run_onecase() {
  local inter_sched=$1
  local run_type=$2
  local seed=$3

  # Determine interslice scheduling objective name
  local inter_sched_name=""
  if [ $inter_sched -eq 9 ]; then
    inter_sched_name="max_throughput"
  elif [ $inter_sched -eq 91 ]; then
    inter_sched_name="proportional_fairness"
  elif [ $inter_sched -eq 92 ]; then
    inter_sched_name="m_lwdf"
  else
    echo "Invalid argument: inter_sched must be 9/91/92"
    return 1
  fi

  # Determine backlogged or customized
  local speed=0
  if [ $run_type == "customize" ]; then
    speed=40
  elif [ $run_type == "backlogged" ]; then
    speed=12
  else
    echo "Invalid argument: run_type must be 'customize'/'backlogged'"
    return 1
  fi

  ../../LTE-Sim SingleCellWithI 1 ${inter_sched} 1 30 ${seed} ${speed} ./configs/config-${run_type}.json 2> ./results/${run_type}_${inter_sched_name}_${seed}.log > /dev/null &
}

run_experiments() {
  local run_type=$1

  for seed in $(seq 0 2); do
    cp ${TRACE_PATH}/mapping${seed}.config ${TRACE_PATH}/mapping.config
    run_onecase 9 ${run_type} ${seed}
    run_onecase 91 ${run_type} ${seed}
    run_onecase 92 ${run_type} ${seed}
    sleep 1
  done
}

if [ ! -d "./results" ]; then
  mkdir ./results
fi

run_experiments customize
run_experiments backlogged
