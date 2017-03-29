#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim_dev
date = 032917
j=0
for ttx in 0 1
do
    export ttx
    for i in -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.150, 0.175, 0.2, 0.225, 0.25
    do
        export j
        export i
        qsub -N job_"$date"_place_cell_nap_i_inj"$i"_ttx"$ttx"_"$j" -l d_rt=36000 -b y -cwd -V "python simulate_place_cell_with_soma_offset.py $i $ttx $j"
        j=$(($j+1))
    done
done