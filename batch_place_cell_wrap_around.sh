#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim_dev
d="032917"
export d
j=0
for ttx in 0 1
do
    export ttx
    for i in -0.1 -0.075 -0.05 -0.025 0. 0.025 0.05 0.075 0.1 0.125 0.150
    do
        for k in {1..10}
        do
            export j
            export i
            qsub -N job_"$d"_place_cell_nap_i_inj"$i"_ttx"$ttx"_"$j" -l d_rt=36000 -b y -cwd -V "python simulate_place_cell_with_soma_offset.py $i $ttx $j"
            j=$(($j+1))
        done
    done
done