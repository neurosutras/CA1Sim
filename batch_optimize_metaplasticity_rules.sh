#!/bin/bash
cd $HOME/CA1Sim_dev
qsub -pe batch 2 -N optimize_012517_cell01 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 1'
qsub -pe batch 2 -N optimize_012517_cell02 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 2'
qsub -pe batch 2 -N optimize_012517_cell03 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 3'
qsub -pe batch 2 -N optimize_012517_cell04 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 4'
qsub -pe batch 2 -N optimize_012517_cell05 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 5'
qsub -pe batch 2 -N optimize_012517_cell06 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 6'
qsub -pe batch 2 -N optimize_012517_cell07 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 7'
qsub -pe batch 2 -N optimize_012517_cell08 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 8'
qsub -pe batch 2 -N optimize_012517_cell10 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 10'
qsub -pe batch 2 -N optimize_012517_cell11 -l d_rt=36000 -b y -cwd -V 'python optimize_metaplasticity_rule.py 11'
