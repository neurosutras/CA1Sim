#!/bin/bash
cd $HOME/CA1Sim_dev
for i in 1 2 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23
do
    s = 'python optimize_discrete_long_plasticity_rule_032017.py '$i
    qsub -pe batch 4 -N export_discrete_long_032617_cell$i -l d_rt=1200 -b y -cwd -V $s
done
for i in 1 2 3 4 5 6 7
do
    s = 'python optimize_discrete_long_plasticity_rule_spont_032017.py '$i
    qsub -pe batch 4 -N export_discrete_long_032617_cell$i -l d_rt=1200 -b y -cwd -V $s
done