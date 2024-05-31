#!/bin/bash 

#SBATCH --job-name=gamlss_slope
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=02:15:00
#SBATCH --account=p274
#SBATCH --output scripts/reliability_scripts/logs/slurm-%j.txt
##SBATCH --partition=hugemem


######################
# setting environment
######################
echo "SETTING UP COLOSSUS ENVIRONMENT"
echo "LOADING SINGULARITY MODULE"
module purge
module load R/4.1.0-foss-2021a
echo `which R`
#module load matlab
#echo `which matlab`

#echo "SOURCING FREESURFER"
#export FREESURFER_HOME=/cluster/projects/p274/tools/mri/freesurfer/current
#source $FREESURFER_HOME/SetUpFreeSurfer.sh
#echo "SOURCING FSL"
#FSLDIR=/cluster/projects/p274/tools/mri/fsl/current
#. ${FSLDIR}/etc/fslconf/fsl.sh
#PATH=${FSLDIR}/bin:${PATH}
#export FSLDIR PATH
export LANG=en_US.utf8


outdir=$1
outrel=$2


echo "$outdir $model $phenotype"

if [ -e scripts/reliability_scripts/logs/slurm.gamlss.slope.txt ]; then
	exit 1
fi

mv scripts/reliability_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/reliability_scripts/logs/slurm.gamlss.slope.txt


Rscript scripts/reliability_scripts/gamlss_slope.r $outdir $outrel
