#!/bin/bash
# 
#SBATCH --job-name="HWinvest"
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:59:00
#SBATCH --partition=esm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=n.wagner@fz-juelich.de
#SBATCH --account=esmtst
#SBATCH --output=HWinvest-out.%j
#SBATCH --error=HWinvest-err.%j


source /p/project/cslts/local/juwels/env_ini.JUWELS.stage2020.GCC
cd /p/scratch/cjibg35/tsmpforecast/development/postpro_genericVAlidationTool/examples

srun -n 48 python examples_DetectHeatwaves_Domain.py

wait

