#!/bin/sh
# You must specify a valid email address!
#SBATCH --mail-user=noelle.schenk@students.unibe.ch
# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="angsd_gbs"

# Runtime and memory
#SBATCH --mem-per-cpu=20G
#SBATCH --time=04:00:00

# For parallel jobs
#SBATCH --cpus-per-task=12

#SBATCH --output=new_restricted01_%j.out
#SBATCH --error=new_restricted01_%j.err

############################# execute code here################################
#NOTE : change error/outfile names, name of -out in angsd line and parameters to change.

# load the required modules
module load vital-it
module load UHTS/Analysis/ANGSD/0.911

angsd -doCounts 1 -GL 2 -minInd 50 -out new_restricted01 -nThreads 12 -doGlf 2 -postCutoff 0.2 -doMajorMinor 1 -doMaf 2 -doPost 1 -doGeno 3 -geno_minDepth 4 -SNP_pval 0.10000 -minMaf -1.000000 -rmTriallelic 0.000000 -bam bamtest.filelist

# -minInd 100 : minimum 100 of the 201 individuals 
# -postCutoff 0.333333 :  (Only genotype to missing if below this threshold)
# -geno_minDepth 5 : -1 is no cutoff
# -minMaf 0.0 : perform a SNP calling, only take variable sites. Remove sites with MAF below ; only print those sites with an allele
#    frequency above 0.05
# todo : either minmaf 0 OR snp_pval 1
# -SNP_pval 1 : remove sites with pvalue > x

# "fixed parameters"
# -doCounts 1 : count stuff, required for other stuff
# -doMajorMinor 1 : input sequencing data like bam files, major and minor allele can be inferred directly from likelihoods with maximum likelihood method
    # required to calculate -doMaf from genotype likelihoods.
# doMaf 2 : 2 is for known major, unknown minor. 1 is for both known (inferred from GL) and 8 is estimated based on base counts.
# -GL 2 : use GATK
# rmTriallelic : remove triallelic situations
# -doGlf 2 : output geno file
# -doMaf 4 (?) Frequency from genotype probabilities
# -doPost 1 : calculate the posterior probability of genotypes with frequency as prior
# -doGeno 11 : call genotypes, 1 : write major and minor ; 2: write the called genotype encoded as -1,0,1,2, -1=not called
    # 8: write the posterior probability of all possible genotypes

# http://www.jamesandthegiantcorn.com/2017/03/26/correcting-genotyping-errors-when-constructing-genetic-maps-from-genotyping-by-sequencing-gbs-data/
# https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
 