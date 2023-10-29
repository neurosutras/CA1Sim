#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export JOB_NAME=CA1Sim_DC_offset_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.o
#SBATCH -e /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.e
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -t 24:00:00
#SBATCH --mail-user=milstein@cabm.rutgers.edu
#SBATCH --mail-type=ALL

set -x

cd $WORK/CA1Sim

ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py 0 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py 1 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py 1 0.2 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset.py 1 0.22 $SCRATCH/data/CA1Sim
wait
EOT