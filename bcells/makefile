%/countfiles : %.tuples.5.1.counts %.tuples.6.2.counts %.tuples.7.2.counts %.tuples.7.3.counts %.tuples.8.2.counts %.tuples.8.3.counts %.tuples.9.3.counts

%/tuples.9.3.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.8.3.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.8.2.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.7.3.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.7.2.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.6.2.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.5.1.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 1 -r 1 -i $< > $@

%/tuples.5.2.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 5 -l 2 -r 2 -i $< > $@

%/tuples.3.1.counts : %.fasta
	python ../tuple-counting/count-paired-tuples.py -w 3 -l 1 -r 1 -i $< > $@


