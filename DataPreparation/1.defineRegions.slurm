## load modules
module load tabix/0.2.6-fasrc01


## create dirs
mkdir 2.1.call
mkdir 2.2.filter
mkdir 2.3.merge


mkdir 2.1.call/heradem
mkdir 2.2.filter/heradem
mkdir 2.3.merge/heradem




#################### Produce Datasets #########################
## variables
# refnam="hmelv25"
# fasref="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/hmelv25/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa"
# mindep=20
refnam="heradem"
fasref="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/heradem/Heliconius_erato_demophoon_v1_-_scaffolds.fa"
mindep=12

export refnam=$refnam
export fasref=$fasref
export mindep=$mindep

mkdir 2.1.call/${refnam}
mkdir 2.2.filter/${refnam}
mkdir 2.3.merge/${refnam}

## 2.1. GT call (copy necessary files)
# for chrom in `cat /n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/$refnam/$refnam.scaff_in_chrom.txt`; do 
#     echo $chrom;
#     for ind in `cat /n/holyscratch01/mallet_lab/fseixas/4.DeepPhylogeny/sppList.list`; do
#         echo $ind;
#         rsync -av --progress /n/mallet_lab/fseixas/1.projects/3.ziheng_yang/0.data/2.1.call/$refnam/$chrom.$ind.to$refnam.vcf.gz 2.1.call/$refnam/
#     done
# done

## 2.2. filter individuals
cat support/mindp.$refnam.minDP$mindep.txt | egrep SAR | xargs -n 2 sh -c '
    sbatch /n/home12/fseixas/code/heliconius_seixas/5.deepPhylogeny/2.2.filter.slurm $0 $refnam $fasref $1 $mindep
'

## 2.3. merge & simplify
sbatch ~/code/heliconius_seixas/5.deepPhylogeny/2.3.merge.slurm $refnam $mindep



#################### Create Alignments #########################
# variables
refnam="hmelv25"
repeat="/n/mallet_lab/fseixas/1.projects/0.basic_data/reference_genomes/${refnam}/${refnam}_repeats.txt"
mindep=12
setloc="exonic.minL100.maxLnan.minGnan"

# 
cd 3.datasets
# 3.0. regions to samples
bash /n/home12/fseixas/code/heliconius_seixas/5.deepPhylogeny/3.0.selectRegions.sh $refnam 
# 3.1. create alignments
sbatch /n/home12/fseixas/code/heliconius_seixas/5.deepPhylogeny/3.1.produceAlignments.slurm $refnam $setloc $repeat $mindep


#ls hmelv25/minDP20/exonic.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP20/exonic.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP20.exonic.loci.txt' 
#ls hmelv25/minDP20/noncod.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP20/noncod.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP20.noncod.loci.txt' 
#ls hmelv25/minDP12/exonic.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP12/exonic.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP12.exonic.loci.txt' 
#ls hmelv25/minDP12/noncod.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP12/noncod.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP12.noncod.loci.txt' 
#ls hmelv25/minDP8/exonic.*/2.vcf2Fasta/  | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP8/exonic.*/2.vcf2Fasta/$0.$1.$2*phy  | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP08.exonic.loci.txt' 
#ls hmelv25/minDP8/noncod.*/2.vcf2Fasta/  | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat hmelv25/minDP8/noncod.*/2.vcf2Fasta/$0.$1.$2*phy  | sed "'"s,ERAtog,ERA,g"'" >> hmelv25.minDP08.noncod.loci.txt' 

#ls heradem/minDP20/exonic.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP20/exonic.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP20.exonic.loci.txt' 
#ls heradem/minDP20/noncod.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP20/noncod.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP20.noncod.loci.txt' 
#ls heradem/minDP12/exonic.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP12/exonic.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP12.exonic.loci.txt' 
#ls heradem/minDP12/noncod.*/2.vcf2Fasta/ | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP12/noncod.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP12.noncod.loci.txt' 
#ls heradem/minDP8/exonic.*/2.vcf2Fasta/  | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP8/exonic.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP08.exonic.loci.txt' 
#ls heradem/minDP8/noncod.*/2.vcf2Fasta/  | grep phy | awk -F'.' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | xargs -n 3 sh -c 'cat heradem/minDP8/noncod.*/2.vcf2Fasta/$0.$1.$2*phy | sed "'"s,ERAtog,ERA,g"'" >> heradem.minDP08.noncod.loci.txt' 


