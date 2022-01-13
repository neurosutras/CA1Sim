#!/bin/bash -l

export j=0
for i in -0.6 -0.4 -0.2 0 0.2 0.4 0.6
do
  export i
  python3 20220112_simulate_place_cell_subtr_inh_record_spines_DC_offset.py 0 $j $i $DATA_DIR &
  ((j=j+1))
  export j
  python3 20220112_simulate_place_cell_subtr_inh_record_spines_DC_offset.py 1 $j $i $DATA_DIR &
  ((j=j+1))
  export j
done
wait