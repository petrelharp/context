simfiles = $(shell ls ising-sims/*.RData | sed -e "s/RData/Rout/g")

.PHONY : analyses

analyses : $(simfiles)

%.Rout : %.RData
	Rscript ising-inference.R -i $< -n 1000 > $*.log 2> $@

test : check-ising.R
	Rscript --vanilla $<

info : 
	Rscript ../siminfo.R ising-sims
