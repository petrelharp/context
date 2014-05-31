
.PHONY : mm9rn5.data mm9rn5

11_filter_regions.sorted.pos : 10_filter_regions.pos
	tail -n +2 10_filter_regions.pos | cut -f 1 | sed -e 's/\(.*\):\([0-9]*\)-\([0-9]*\).*/\1 \2 \3/' | sort > $@

mm9rn5counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9rn5/*))

mm9rn5.data :
	-mkdir mm9rn5
	cd mm9rn5; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsRn5/axtNet/chr$${x}.mm9.rn5.net.axt.gz; done
	cd mm9rn5; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsRn5/axtNet/chrX.mm9.rn5.net.axt.gz

mm9rn5 : $(mm9rn5counts)
	echo $(mm9rn5counts)

%.axt.sub.gz : %.axt.gz
	python ../tuple-counting/get-regions-from-alignment.py --axtfile $< --posfile 11_filter_regions.sorted.pos -o $@

%.5.1.counts : %.net.axt.sub.gz
	python ../tuple-counting/count-paired-tuples.py --strict --infile $< --informat axt -w 5 -l 1 -r 1 --outfile $@
