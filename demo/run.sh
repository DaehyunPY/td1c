#!/bin/bash
#PBS -l walltime=2400:00:00
#PBS -l nodes=1:ppn=12

export OMP_NUM_THREADS=12
export OMP_SCHEDULE=static
export OMP_MAX_ACTIVE_LEVELS=1

IFS=$'\n'

dir="${PBS_O_WORKDIR:-$(cd "$(dirname "$0")"; pwd -P)}"
inputs=($(find "$dir" -name "*.inp" -maxdepth 1 2> /dev/null))
for inp in "${inputs[@]}"
do
    log="${inp%.inp}.log"
    [[ -f "$inp" ]] && [[ ! -e "$log" ]] && td1c "$inp" &> "$log"
done
