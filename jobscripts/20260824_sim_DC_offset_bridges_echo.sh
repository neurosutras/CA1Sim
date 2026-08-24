#!/bin/bash -l

source $HOME/cpu_py311_intelmpi.sh

cd $PROJECT/CA1Sim

declare seed=0

for ((i=-2; i<8; i++))
do
  current=$(awk -v i="$i" 'BEGIN { printf "%.3f", i * 0.026 }')
  for ((j=0; j<10; j++))
  do
    echo -n 1 python 20231029_simulate_place_cell_record_syn_currents_DC_offset.py $seed $current $SCRATCH/data/CA1Sim &
    ((++seed))
  done
done
wait
