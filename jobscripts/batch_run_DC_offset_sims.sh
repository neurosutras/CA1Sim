#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export LABEL=CA1Sim_DC_offset
export JOB_NAME="$LABEL"_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.o
#SBATCH -e /scratch1/06441/aaronmil/logs/CA1Sim/$JOB_NAME.%j.e
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -t 1:00:00
#SBATCH --mail-user=neurosutras@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -x

cd $WORK2/CA1Sim

export DATA_DIR=$SCRATCH/data/CA1Sim

declare -a offset=(-0.6 -0.4 -0.2 0 0.2 0.4 0.6)

arraylength=${#offset[@]}

for ((i=0; i<${arraylength}; i++))
do
  ibrun -n 1 python3 20220112_simulate_place_cell_subtr_inh_record_spines_DC_offset.py 0 $i ${offset[$i]} $DATA_DIR &
  let "j = $i + $arraylength"
  ibrun -n 1 python3 20220112_simulate_place_cell_subtr_inh_record_spines_DC_offset.py 1 $j ${offset[$i]} $DATA_DIR &
done
EOT
