%.Rout : %.RData
	Rscript --vanilla ising-inference.R -i $< > $*.log 2> $@

test : check-ising.R
	Rscript --vanilla $<
