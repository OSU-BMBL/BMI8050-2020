#!/usr/bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=16:gpus=1
#PBS -l mem=24GB
#PBS -m n
#PBS -j oe
#PBS -S /usr/bin/bash
#PBS -A PAS1791

module load python/3.6 cuda/9.2.88
cd /fs/ess/PAS1791/course_code
echo 'Hello world'
