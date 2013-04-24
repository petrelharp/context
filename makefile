.PHONY : test

chromsomes = 2L 2R 3L 3R 4 X
codons = GCT GCC GCA GCG CGT CGC CGA CGG AGA AGG AAT AAC GAT GAC TGT TGC CAA CAG GAA GAG GGT GGC GGA GGG CAT CAC ATT ATC ATA TTA TTG CTT CTC CTA CTG AAA AAG ATG TTT TTC CCT CCC CCA CCG TCT TCC TCA TCG AGT AGC ACT ACC ACA ACG TGG TAT TAC GTT GTC GTA GTG
types = orfcoding nonorfcoding noncoding

dmel-all-r5.50.gff.gz :
	wget ftp://ftp.flybase.net/releases/current/dmel_r5.50/gff/dmel-all-r5.50.gff.gz

dmel_mitochondrion_genome.raw.gz YHet.raw.gz XHet.raw.gz X.raw.gz Uextra.raw.gz U.raw.gz 4.raw.gz 3RHet.raw.gz 3R.raw.gz 3LHet.raw.gz 3L.raw.gz 2RHet.raw.gz 2R.raw.gz 2LHet.raw.gz 2L.raw.gz : 
	wget -R ftp://ftp.flybase.net/releases/current/dmel_r5.50/dna/

dsim-all-chromosome-r1.4.fasta.gz :
	wget ftp://ftp.flybase.net/releases/current/dsim_r1.4/fasta/dsim-all-chromosome-r1.4.fasta.gz

dyak-all-chromosome-r1.3.fasta.gz :
	wget ftp://ftp.flybase.net/releases/current/dyak_r1.3/fasta/dyak-all-chromosome-r1.3.fasta.gz

dmel-all-chromosome-r5.50.fasta.gz :
	wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.50_FB2013_02/fasta/dmel-all-chromosome-r5.50.fasta.gz

# $(chromosomes).raw.$(codons).gz : get-all-positions.sh get-positions.py
# 	source @<

CDS-dmel-$(chromosomes)-r5.50.CDS.starts.ends.gz : extract-coding.sh
	source @<

# $(chromosomes).raw.$(codons).$(orfcoding).gz : is-coding.sh is-coding.py
# 	source @<

freeze2.vcf.gz : 
	wget http://www.hgsc.bcm.tmc.edu/projects/dgrp/freeze2_Feb_2013/vcf_files/freeze2.vcf.gz

freeze2.snps.vcf.gz : freeze2.vcf.gz
	# grep for SNPs only
	( zcat freeze2.vcf.gz | head -n 19; zcat freeze2.vcf.gz | grep "^.*\s.*\s.*\s[ACGT]\s[ACGT]\s" ) | gzip -c > freeze2.snps.vcf.gz

freeze2.snps-freqs.vcf.gz : freeze2.snps.vcf.gz 
	# then, hist these and pick cutoffs for next step...
	# add total coverage as first column
	#   and then only keep snps with coverage at least 3000 and no more than 10000 
	# OLD:  ( zcat freeze2.snps.vcf.gz | head -n 19; zcat freeze2.snps.vcf.gz | tail -n +20 | paste coverages - | awk 'int($$1)>=3000 && int($$1)<=10000' | cut -f 2- ) | gzip -c > freeze2.snps.filtered.vcf.gz
	( echo "coverage" ; zcat freeze2.snps.vcf.gz | tail -n +20 | cut -f 10- | sed -e "s/[0-9\/.]*:\([0-9]*\):\([0-9]*\):\([0-9.]*\)/\1 \2 /g" | awk '{ for(i=1; i<=NF;i++) j+=$$i; print j; j=0 }' ) > coverages
	( echo "quality" ; zcat freeze2.snps.vcf.gz | tail -n +20 | cut -f 10- | sed -e "s/[0-9\/.]*:\([0-9]*\):\([0-9]*\):\([0-9.]*\)/\3 /g" | awk '{ for(i=1; i<=NF;i++) j+=$$i; print j; j=0 }' ) > qualities
	pseq freeze2.snps.vcf.gz counts | paste - coverages  | paste - qualities | gzip -c > freeze2.snps-freqs.vcf.gz

freeze2.snps-freqs-meta.vcf.gz : freeze2.snps-freqs.vcf.gz
	pseq freeze2.snps.vcf.gz meta-matrix | cut -f 3- >meta-matrix.temp
	zcat freeze2.snps-freqs.vcf.gz | paste - meta-matrix.temp | gzip -c > freeze2.snps-freqs-meta.vcf.gz

dm3.droYak2.all.chain.gz :
	wget http://hgdownload.cse.ucsc.edu/goldenPath/dm3/vsDroYak2/dm3.droYak2.all.chain.gz

dm3.droSim1.all.chain.gz :
	wget http://hgdownload.cse.ucsc.edu/goldenPath/dm3/vsDroSim1/dm3.droSim1.all.chain.gz

freeze2.snps-dsim-mapped freeze2.snps-dyak-mapped : dm3.droYak2.all.chain.gz dm3.droSim1.all.chain.gz
	zcat  freeze2.snps-freqs.vcf.gz | tail -n +2 | cut -f 1 | sed -e "s/:\(.*\)/:\1-\1/" | gzip -c > freeze2.snps-positions.gz
	../liftOver -positions freeze2.snps-positions.gz ../flybase/dm3.droSim1.all.chain.gz freeze2.snps-dsim-mapped freeze2.snps-dsim-unmapped
	../liftOver -positions freeze2.snps-positions.gz ../flybase/dm3.droYak2.all.chain.gz freeze2.snps-dyak-mapped freeze2.snps-dyak-unmapped


dsim.fasta.gz dyak.fast.gz : 
	zcat dsim-all-chromosome-r1.4.fasta.gz | awk '/^>/ { if (line) print line; print $0 } !/^>/ { line = line $0 } END { print line }' | fold -w 210 | gzip -c > dsim.fasta.gz
	zcat dyak-all-chromosome-r1.3.fasta.gz | awk '/^>/ { if (line) print line; print $0 } !/^>/ { line = line $0 } END { print line }' | fold -w 210 | gzip -c > dyak.fasta.gz

###
# Stuff for checking
#
# Make a list of data frame, the k-th of which has the k-th component in the VCF entries:

test :
	echo "y <- read.table(pipe('zcat freeze2.snps.vcf.gz | head -n 100'),header=FALSE); \
		y <- y[,10:ncol(y)]; \
		options(warn=-1);\
		zzz <- lapply( 2:4, function (k) { data.frame( lapply( y, function (z) { as.numeric( sapply( strsplit(levels(z),':'), '[', k ) )[as.numeric(z)] } ) ) } ); \
		freqs <- read.table(pipe('zcat freeze2.snps-freqs.vcf.gz | head -n 82'), header=TRUE); \
		if (!all( freqs['quality'] == rowSums(zzz[[3]],na.rm=TRUE) ) | !all( freqs['coverage'] == rowSums(zzz[[1]]) + rowSums(zzz[[2]]) ) ) { stop('Oops!') } else { print('Checks out!  All good!') } " | Rscript -

