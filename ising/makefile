%.Rout : %.RData
	Rscript --vanilla ising-inference.R $< > $*.log 2> $@

test : check-ising.R
	Rscript --vanilla $<
