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
	wget ftp://ftp.flybase.net//genomes/Drosophila_melanogaster/dmel_r5.50_FB2013_02/fasta/dmel-all-chromosome-r5.50.fasta.gz

$(chromosomes).raw.$(codons).gz : get-all-positions.sh get-positions.py
	source @<

CDS-dmel-$(chromosomes)-r5.50.CDS.starts.ends.gz : extract-coding.sh
	source @<

$(chromosomes).raw.$(codons).$(orfcoding).gz : is-coding.sh is-coding.py
	source @<

freeze2.vcf.gz : 
	wget http://www.hgsc.bcm.tmc.edu/projects/dgrp/freeze2_Feb_2013/vcf_files/freeze2.vcf.gz

freeze2.snps.vcf.gz freeze2.snps-freqs.vcf.gz : freeze2.vcf.gz coverages.gz
	# grep for SNPs only
	( zcat freeze2.vcf.gz | head -n 19; zcat freeze2.vcf.gz | grep "^.*\s.*\s.*\s[ACGT]\s[ACGT]\s" ) | gzip -c > freeze2.snps.vcf.gz
	# add total coverage as first column
	gunzip coverages.gz
	# only keep snps with coverage at least 3000 and no more than 10000 (see hist of coverages)
	( zcat freeze2.snps.vcf.gz | head -n 19; zcat freeze2.snps.vcf.gz | paste coverages - | awk 'int($1)>=3000 && int($1)<=10000' | cut -f 2- ) | gzip -c > freeze2.snps.filtered.vcf.gz
	gzip coverages
	pseq freeze2.snps.filtered.vcf.gz counts | gzip -c > freeze2.snps-freqs.vcf.gz

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

coverages.gz : freeze2.snps.vcf.gz
	zcat freeze2.snps.vcf.gz | tail -n +20 | cut -f 11- | sed -e "s/[0-9\/.]*:\([0-9]*\):\([0-9]*\):[0-9.]*/\1 \2 /g" | awk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' | gzip -c > coverages.gz
