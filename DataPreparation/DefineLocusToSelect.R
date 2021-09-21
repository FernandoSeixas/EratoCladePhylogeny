## load libraries ----
require(optparse)

option_list <- list(
		make_option(c('-d','--direc'), action='store', type='character', default=NULL, help='path to directory containing input file'),
		make_option(c('-i','--input'), action='store', type='character', default=NULL, help='input file name'),
		make_option(c('-o','--output'),action='store', type='character', default=NULL, help='output file name'),
		make_option(c('-l','--minl'),	action='store', type='numeric', default=100,   help='minimum locus size'),
		make_option(c('-m','--maxl'),	action='store', type='numeric', default=1000,  help='maximum locus size'),
		make_option(c('-g','--gap'),	action='store', type='numeric', default=2000,  help='minimum distance between loci')
)

opt <- parse_args(OptionParser(option_list = option_list))

## read files and transform
dir <- opt$direc
inp <- read.table(paste0(dir, opt$input))
names(inp) = c("scaffold","sta","end")
inp$sta = ifelse(inp$sta == 0, 1, inp$sta)
inp$len = inp$end - inp$sta + 1

## define variables
minlen = opt$minl
maxlen = opt$maxl
mingap = opt$gap


## subset intervals by size
fltinp = inp
fltinp = subset(fltinp, len >= minlen)
print(nrow(fltinp))

## define functions
splitseg <- function(x) {
	spl = data.frame()
	newend = x$sta + maxlen - 1
	df = data.frame(scaffold = x$scaffold, sta = x$sta, end = newend)
	spl = rbind(spl, df)
	while (newend < x$end) {
		newsta = newend + mingap + 1
		newend = newsta + maxlen - 1
		if (newend > x$end) {newend = x$end}
		if (newsta < x$end & (x$end-newsta) > minlen) {
			df = data.frame(scaffold = x$scaffold, sta = newsta, end = newend)
			spl = rbind(spl, df)
		} 
	}
	spl$len = spl$end - spl$sta + 1
	return(spl)
}


## subset regions according to filters [minlen, maxlen, gapslen]
regions_set = data.frame()

for (chrom in 1:length(unique(fltinp$scaffold))) {
	print(unique(fltinp$scaffold)[chrom])
	sub = subset(fltinp, scaffold == unique(fltinp$scaffold)[chrom])
	reg = data.frame()
	
	##	for first segment first segment ....
	fir_seg <- sub[1,]
	# if segment within range
	reg = fir_seg
	# if first segment larger than maxlen but can only fit one segment
	if (fir_seg$len > maxlen & fir_seg$len < (2*maxlen+mingap)) {
		fir_seg$end = fir_seg$sta + maxlen - 1
		fir_seg$len = fir_seg$end - fir_seg$sta + 1
		reg = fir_seg
	}
	# split first segment if necessary
	if (fir_seg$len >= (2*maxlen+mingap)) {reg = splitseg(fir_seg)}
	# redefine first segment len (if needed)
	reg$len = reg$end - reg$sta + 1
	
	## for remaining segments ....
	if (nrow(sub) > 1) {
		for (rw in 2:nrow(sub)) {
			seg = sub[rw,]
			ms = reg$end[nrow(reg)]+mingap
			# end of new region more than ##-kb appart and can contain a segment larger than minlen
			if (seg$end > ms & (seg$end-ms) > minlen ) {
				# update sta position if < min start (accomodating gap)
				if (seg$sta < ms) {seg$sta = ms}
				seg$len = seg$end-seg$sta
				# if segment within range
				if (seg$len < maxlen) {
					reg = rbind(reg, seg)
				}
				# if segment > range but cannot accomodate more than two loci
				if (seg$len > maxlen & seg$len < (2*maxlen+mingap)) {
					seg$end = seg$sta + maxlen - 1
					seg$len = seg$end - seg$sta + 1
					reg = rbind(reg, seg)
				}
				# if segment can accomodate more than 1 loci
				if (seg$len >= (2*maxlen+mingap)) {reg = rbind(reg, splitseg(seg))}
			}
		}
	}
	regions_set = rbind(regions_set, reg)
}

write.table(file=paste0(dir, opt$out), regions_set[,1:3], row.names=F, col.names=F, quote=F, sep="\t")

