%.Rout : %.RData
	Rscript ising-inference.R -i -n 1000 $< > $*.log 2> $@

test : check-ising.R
	Rscript --vanilla $<
