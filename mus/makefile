

.PHONY : mm9rn5.data mm9rn5

11_filter_regions.sorted.pos : 10_filter_regions.pos
	tail -n +2 10_filter_regions.pos | cut -f 1 | sed -e 's/\(.*\):\([0-9]*\)-\([0-9]*\).*/\1 \2 \3/' | sort > $@

%.axt.sub.gz : %.axt.gz
	python ../tuple-counting/get-regions-from-alignment.py --axtfile $< --posfile 11_filter_regions.sorted.pos -o $@

%.5.1.counts : %.net.axt.sub.gz
	python ../tuple-counting/count-paired-tuples.py --strict --infile $< --informat axt -w 5 -l 1 -r 1 --outfile $@

%/makefile :
	@echo '# automatically generated' > $@
	@echo '' >> $@
	@echo 'counts = $$(patsubst %.net.axt.gz,%.5.1.counts,$$(wildcard *.net.axt.gz))' >> $@
	@echo '' >> $@
	@echo '5.1.counts : $$(counts)' >> $@
	@echo '	python ../../tuple-counting/sum-counts.py -o $$@ -i $$<' >> $@
	@echo '' >> $@
	@echo 'include ../makefile' >> $@

############# WET makefile, now redundant

# mouse-rat
mm9rn5.data :
	# mouse-rat
	-mkdir mm9rn5
	cd mm9rn5; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsRn5/axtNet/chr$${x}.mm9.rn5.net.axt.gz; done
	cd mm9rn5; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsRn5/axtNet/chrX.mm9.rn5.net.axt.gz
	make mm9rn5/makefile

mm9rn5counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9rn5/*))

mm9rn5/5.1.counts : $(mm9rn5counts)
	python ../tuple-counting/sum-counts.py -o $@ -i $<

mm9rn5 : $(mm9rn5counts)
	echo $(mm9rn5counts)

# mouse-human
mm9hg19.data :
	# mouse-human
	-mkdir mm9hg19
	cd mm9hg19; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsHg19/axtNet/chr$${x}.mm9.hg19.net.axt.gz; done
	cd mm9hg19; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsHg19/axtNet/chrX.mm9.hg19.net.axt.gz

mm9hg19counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9hg19/*))

mm9hg19/5.1.counts : $(mm9hg19counts)
	python ../tuple-counting/sum-counts.py -o $@ -i $<

mm9hg19 : $(mm9hg19counts)
	echo $(mm9hg19counts)

# mouse-platypus
mm9ornAna1counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9ornAna1/*))

mm9ornAna1.data :
	# mouse-platypus
	-mkdir mm9ornAna1
	cd mm9ornAna1; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsOrnAna1/axtNet/chr$${x}.mm9.ornAna1.net.axt.gz; done
	cd mm9ornAna1; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsOrnAna1/axtNet/chrX.mm9.ornAna1.net.axt.gz

mm9ornAna1 : $(mm9ornAna1counts)
	echo $(mm9ornAna1counts)

# mouse-guinea pig
mm9oryCun2counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9oryCun2/*))

mm9oryCun2.data :
	# mouse-guinea pig
	-mkdir mm9oryCun2
	cd mm9oryCun2; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsOryCun2/axtNet/chr$${x}.mm9.oryCun2.net.axt.gz; done
	cd mm9oryCun2; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsOryCun2/axtNet/chrX.mm9.oryCun2.net.axt.gz

mm9oryCun2 : $(mm9oryCun2counts)
	echo $(mm9oryCun2counts)

# mouse-chicken
mm9galGal3counts = $(patsubst %.net.axt.gz,%.5.1.counts,$(wildcard mm9galGal3/*))

mm9galGal3.data :
	# mouse-chicken
	-mkdir mm9galGal3
	cd mm9galGal3; for x in $$(seq 19); do wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsGalGal3/axtNet/chr$${x}.mm9.galGal3.net.axt.gz; done
	cd mm9galGal3; wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/vsGalGal3/axtNet/chrX.mm9.galGal3.net.axt.gz

mm9galGal3 : $(mm9galGal3counts)
	echo $(mm9galGal3counts)
