#!/bin/bash
#SBATCH -n 8    # Number of cores requested
#SBATCH -N 1    # Ensure that all cores are on one machine
#SBATCH -t 60 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o prepGPhocs_%j.out # Standard out goes to this file
#SBATCH -e prepGPhocs_%j.err # Standard err goes to this filehostname


## load modules
module load bedtools2/2.26.0-fasrc01
module load tabix/0.2.6-fasrc01
module load R/3.6.3-fasrc01
export R_LIBS_USER=$HOME/apps/R_v3.6.3:$R_LIBS_USER


## prepare global dataset **********************************************************************************************

## variables & export
refnam="heradem"; export refnam=$refnam; 
minlen=100;       export minlen=$minlen;
maxlen=250;       export maxlen=$maxlen;
gaplen=2000;      export gaplen=$gaplen;
buffer=2000;      export buffer=$buffer;

## Define path to scaffold sizes and genome annotations
if [ $refnam = "hmelv25" ];
then
    chromoSizes="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/hmelv25.scaffold.sizes";
    gff="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/Heliconius_melpomene_melpomene_Hmel2.5.gff3"
fi
if [ $refnam = "heradem" ];
then
    chromoSizes="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/heradem.scaffold.sizes";
    gff="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/Heliconius_erato_demophoon_v1.gff3"
fi


## Exon coordinates file
grep "CDS" ${gff} | cut -f1,4,5 | sed 's, ,\t,g' | sort -k1,1 -k2,2n | uniq > /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/hmelv25.exons_coord.bed
grep "CDS" ${gff} | cut -f1,4,5 | sed 's, ,\t,g' | sort -k1,1 -k2,2n | uniq > /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/heradem.exons_coord.bed

if [ $refnam = "hmelv25" ];
then
    exonicRegions="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/hmelv25.exons_coord.bed";
fi
if [ $refnam = "heradem" ];
then
    exonicRegions="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/heradem.exons_coord.bed";
fi


## Create bed files for NONCOD (intergenic + intronic) and CODING (exonic) regions *************************************

## NONCOD //////
# get non-exonic regions coordinates [complement exonic regions]
bedtools complement -i $exonicRegions -g $chromoSizes > $refnam.noncod.Step1.bed
# run Rscript to get coordinates of valid regions to sample from
Rscript \
    DefineLocusToSelect.R \
    -d ./ \
    -i $refnam.noncod.Step1.bed \
    -l $minlen \
    -m $maxlen \
    -g $gaplen \
    -b $buffer \
    -o $refnam.noncod.Step2.bed

## EXONIC //////
# get exonic unique coordinates (sorted), and only in chromosomes
cat $exonicRegions > $refnam.exonic.Step1.bed
# run Rscript to get coordinates of valid regions to sample from
Rscript \
    DefineLocusToSelect.R \
    -d ./ \
    -i $refnam.exonic.Step1.bed \
    -l $minlen \
    -m $maxlen \
    -g $gaplen \
    -b $buffer \
    $refnam.exonic.Step2.bed