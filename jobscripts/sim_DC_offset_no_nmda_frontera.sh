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

ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 20 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 21 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 22 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 23 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 24 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 25 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 26 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 27 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 28 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 29 0.0 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 30 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 31 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 32 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 33 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 34 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 35 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 36 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 37 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 38 0.18 $SCRATCH/data/CA1Sim &
ibrun -n 1 python3 20231029_simulate_place_cell_record_syn_currents_DC_offset_no_nmda.py 39 0.18 $SCRATCH/data/CA1Sim

wait
EOT
