#!/bin/bash

TRACE_PATH="../../cqi-traces-noise0"

run_onecase() {
  local inter_sched=$1
  local inter_sched_name=""

  if [ "$inter_sched" -eq 9 ]; then
    inter_sched_name="max_throughput"
  elif [ "$inter_sched" -eq 91 ]; then
    inter_sched_name="proportional_fairness"
  elif [ "$inter_sched" -eq 92 ]; then
    inter_sched_name="m_lwdf"
  else
    echo "Invalid argument: inter_sched must be 9/91/92"
    return 1
  fi

  ../../LTE-Sim SingleCellWithI 1 "$inter_sched" 1 30 "$seed" 40 ./config.json 2> "./${inter_sched_name}_${seed}.log" > /dev/null &
}

run_experiments() {
  for seed in $(seq 0 2); do
    cp "${TRACE_PATH}/mapping${seed}.config" "${TRACE_PATH}/mapping.config"
    run_onecase 9
    run_onecase 91
    run_onecase 92
    sleep 1
  done
}

run_experiments
