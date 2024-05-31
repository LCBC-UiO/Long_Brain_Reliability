#!/bin/bash 

#SBATCH --job-name=computational_reliability
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=24:00:00
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
deltafile=$2
errorfile=$3
intfile=$4
modelfile=$5
roifile=$6
n_subjects=$7
n_icc=$8

echo $outdir $deltafile $errorfile $intfile $modelfile $roifile $n_subjects $n_icc
oo=$(basename $outdir)
if [ -e scripts/reliability_scripts/logs/slurm.$oo.txt ]; then
	exit 1
fi

mv scripts/reliability_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/reliability_scripts/logs/slurm.$oo.txt


Rscript scripts/reliability_scripts/computational_reliability_aggregate.r $outdir $deltafile $errorfile $intfile $modelfile $roifile $n_subjects $n_icc