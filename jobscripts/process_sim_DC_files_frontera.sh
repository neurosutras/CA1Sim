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
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-user=milstein@cabm.rutgers.edu
#SBATCH --mail-type=ALL

set -x

cd $WORK/CA1Sim

ibrun -n 1 python3 process_and_merge_i_syn_DC_offset_data.py \
  --source-data-config-file-path=config/20231101_i_syn_DC_offset_file_names.yaml \
  --source-data-dir=$SCRATCH/data/CA1Sim \
  --target-file-path=data/20231101_merged_i_syn_DC_offset_data.hdf5
EOT
