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

freeze2.snps.vcf.gz : freeze2.vcf.gz
	( zcat freeze2.vcf.gz | head -n 19; zcat freeze2.vcf.gz | grep "^.*\s.*\s.*\s[ACGT]\s[ACGT]\s" ) | gzip -c > freeze2.snps.vcf.gz
