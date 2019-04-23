#!/bin/sh
# You must specify a valid email address!
#SBATCH --mail-user=noelle.schenk@students.unibe.ch
# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

#SBATCH --job-name="map"

# Mandatory resources (h_cpu=hh:mm:ss)
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00

# Only reserve scratch space if you really need it
#SBATCH --tmp=8G

# Which parallel_environment to use (smp or mpi) and how many slots
#SBATCH --cpus-per-task=4

# Other options

# change
# Array job containing 10 tasks, submitting 4 simultaneously
#SBATCH --array=4-96%4 

# working directory
# %j is the job ID and %a is the array job number (array task ID)
# the output and error files are in the directory where the script is run.
#SBATCH --workdir=.

# change
#SBATCH --output=f7.plate4.%a.%j.out
#SBATCH --error=f7.plate4.%a.%j.err

############################# execute code here################################ 
READDIR='/gpfs/homefs/ips/ns12a545/gbs_to_ref/reads/hq_reads'

# load the required modules
module load vital-it
module load UHTS/Aligner/bowtie2/2.3.0

# indexing the reference genome
# bowtie2-build path-to-genome/genome genome_output_name

bowtie2 -p 4 -x P.EXSERTA.contigs.v1.1.3 -U ${READDIR}/plate4/plate4.${SLURM_ARRAY_TASK_ID}.q30l50.fq -S f7.plate4.${SLURM_ARRAY_TASK_ID}.q30l50.sam

samtools view -bS f7.plate4.${SLURM_ARRAY_TASK_ID}.q30l50.sam | samtools sort - f7.plate4.${}
samtools index

# -x index_name -U reads -S output_sam_name

