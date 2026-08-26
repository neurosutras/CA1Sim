#!/bin/bash -l
#SBATCH -J CA1Sim_DC_offset
#SBATCH -o /expanse/lustre/scratch/aaronmil/temp_project/logs/CA1Sim/CA1Sim_DC_offset.%j.o
#SBATCH -e /expanse/lustre/scratch/aaronmil/temp_project/logs/CA1Sim/CA1Sim_DC_offset.%j.e
#SBATCH -p compute
#SBATCH -N 7
#SBATCH --ntasks-per-node=16
#SBATCH -n 100
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem=249208M
#SBATCH --account=sua199
#SBATCH --export=ALL
#SBATCH --mail-user=milstein@cabm.rutgers.edu
#SBATCH --mail-type=ALL
#SBATCH --constraint="lustre"
#SBATCH --no-requeue

source $HOME/cpu_py311_intelmpi.sh

cd $PROJECT/CA1Sim

MEM_PER_CPU="2492M"

declare seed=0

for ((i=-2; i<8; i++))
do
  current=$(awk -v i="$i" 'BEGIN { printf "%.3f", i * 0.026 }')
  for ((j=0; j<10; j++))
  do
    srun --nodes=1 --ntasks=1 -c 8 --mem-per-cpu=$MEM_PER_CPU --exact --exclusive python \
      20231029_simulate_place_cell_record_syn_currents_DC_offset.py $seed $current $SCRATCH/data/CA1Sim &
    ((++seed))
  done
done
wait
