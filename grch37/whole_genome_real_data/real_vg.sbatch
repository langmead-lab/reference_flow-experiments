#!/bin/bash
#SBATCH --job-name=vg
#SBATCH --time=24:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=512GB
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=end
#SBATCH --mail-user=cnaechy1@jhu.edu

/home-1/cnaechy1@jhu.edu/bin/time-1.9 -v -o /storage2/naechyun/refflow/wg/vg/wg-vg_index_slurm.time_log /net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/software/vg/vg index -g /storage2/naechyun/refflow/wg/vg/wg-0.1.gcsa -k 16 -t 16 -b /storage2/naechyun/refflow/wg/vg/tmp /storage2/naechyun/refflow/wg/vg/wg-0.1.vg -Z 8192
