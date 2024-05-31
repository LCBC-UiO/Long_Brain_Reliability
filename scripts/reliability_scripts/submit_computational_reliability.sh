#!/bin/bash 

#SBATCH --job-name=computational_reliability
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
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
analysis=$2
feature=$3
oname=$4
meanF=$5
seF=$6
meanD=$7
seD=$8
meanE=$9
n_subjects=${10}
n_icc=${11}


echo $outdir $analysis $feature $oname $meanF $seF $meanD $seD $meanE $n_subjects $n_icc

if [ -e scripts/reliability_scripts/logs/slurm.$feature.$oname.txt ]; then
	exit 1
fi

mv scripts/reliability_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/reliability_scripts/logs/slurm.$feature.$oname.txt


Rscript scripts/reliability_scripts/computational_reliability.r $outdir $analysis $feature $oname $meanF $seF $meanD $seD $meanE $n_subjects $n_icc