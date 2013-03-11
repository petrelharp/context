chromsomes = 2L 2R 3L 3R 4 X
codons = GCT GCC GCA GCG CGT CGC CGA CGG AGA AGG AAT AAC GAT GAC TGT TGC CAA CAG GAA GAG GGT GGC GGA GGG CAT CAC ATT ATC ATA TTA TTG CTT CTC CTA CTG AAA AAG ATG TTT TTC CCT CCC CCA CCG TCT TCC TCA TCG AGT AGC ACT ACC ACA ACG TGG TAT TAC GTT GTC GTA GTG
types = orfcoding nonorfcoding noncoding

$(chromosomes).raw.$(codons).gz : get-all-positions.sh get-positions.py
	source @<

CDS-dmel-$(chromosomes)-r5.50.CDS.starts.ends.gz : extract-coding.sh
	source @<

$(chromosomes).raw.$(codons).$(orfcoding).gz : is-coding.sh is-coding.py
	source @<

