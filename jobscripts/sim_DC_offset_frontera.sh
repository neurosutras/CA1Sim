#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export JOB_NAME=CA1Sim_DC_offset_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.o
#SBATCH -e /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.e
#SBATCH -p small
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 4:00:00
#SBATCH --mail-user=milstein@cabm.rutgers.edu
#SBATCH --mail-type=ALL

set -x

cd $WORK/CA1Sim

export j=0
for i in {0..9}
do
  export i
  ((j=i+10))
  export j
  ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py $i 0.0 $SCRATCH/data/CA1Sim &
  ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py $j 0.18 $SCRATCH/data/CA1Sim &
done
wait
EOT
