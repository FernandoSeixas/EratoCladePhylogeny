#!/bin/bash
#SBATCH -n 16   # Number of cores requested
#SBATCH -N 1   # Ensure that all cores are on one machine
#SBATCH -t 600 # Runtime in minutes
#SBATCH -p general,shared,serial_requeue # Partition to submit to
#SBATCH --mem-per-cpu 4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o bppaligns_%j.out # Standard out goes to this file
#SBATCH -e bppaligns_%j.err # Standard err goes to this filehostname


# load modules
module load bcftools/1.5-fasrc02


## variables ***********************************************************************************************************
refnam=$1  ## name of the reference genome
claden=$2  ## name of the clade to analyze
setloc=$3  ## name of the set of loci to analyze
mdipDP=$4  ## filter used for the minimum depth of coverage 
repeats=$5 ## file containing the annotated repeats in the reference genomes

export refnam=$refnam
export claden=$claden
export setloc=$setloc
export mdipDP=$mdipDP
export repeats=$repeats

# create directories
mkdir $claden.$mdipDP
cd $claden.$mdipDP
mkdir $setloc
mkdir $setloc/1.vcf_regions
mkdir $setloc/2.phy_iupac

# extract regions (loci) from vcf
cat $refnam.$seloc.Step2.bed  | xargs -n 3 -P 16 sh -c '/n/helmod/apps/centos7/Core/bcftools/1.5-fasrc02/bin/bcftools view -r $0:$1-$2 --targets-file ^$repeats /n/scratchlfs/mallet_lab/fseixas/ziheng_yang/0.data/2.4.norm/$refnam/$claden/$0.mapto$refnam.$claden.$mdipDP.vcf.gz > $setloc/1.vcf_regions/$0_$1_$2.vcf'

# transform to phy format [input .vcf, sta, end, gap/site, gap/loci, output .phy, output .bed]
rm $setloc/$refnam.$setloc.LociPassingFilters.bed
cat $refnam.$seloc.Step2.bed | xargs -n 3 -P 16 sh -c '
    perl vcf2phy.iupac.pl \
    $setloc/1.vcf_regions/$0_$1_$2.vcf $0 $1 $2 0 0.5 \
    $setloc/2.phy_iupac/$0_$1_$2.phy \
    $setloc/$refnam.$setloc.LociPassingFilters.bed
' > $setloc/vcf2phy.log

# concatenate all files
rm $setloc/$refnam.Z.$setloc.loci.txt
rm $setloc/$refnam.A.$setloc.loci.txt
rm $setloc/$refnam.U.$setloc.loci.txt
for chrom in `cat ../$refnam.sexzz.txt`; do for file in `ls $setloc/2.phy_iupac/$chrom*.phy`; do cat $file >> $setloc/$refnam.Z.$setloc.loci.txt; done; done; ## Z-chromosome
for chrom in `cat ../$refnam.autos.txt`; do for file in `ls $setloc/2.phy_iupac/$chrom*.phy`; do cat $file >> $setloc/$refnam.A.$setloc.loci.txt; done; done; ## autosomes
for chrom in `cat ../$refnam.unkno.txt`; do for file in `ls $setloc/2.phy_iupac/$chrom*.phy`; do cat $file >> $setloc/$refnam.U.$setloc.loci.txt; done; done; ## unannotated scaffolds 

# sort loci by coordinates (first scaffold and the start)
sort -k1,1 -k2,2n $setloc/$refnam.$setloc.LociPassingFilters.bed > tmp; mv tmp $setloc/$refnam.$setloc.LociPassingFilters.bed 

# get idea of the size of alignments
ls $setloc/2.phy_iupac/ | xargs -n 1 -P 16 sh -c 'head -n 1 $setloc/2.phy_iupac/$0' > $setloc/$setloc.alignment_sizes.txt
sed -i 's/ /\t/g' $setloc/$setloc.alignment_sizes.txt

