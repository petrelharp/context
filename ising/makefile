simfiles = $(shell ls ising-sims/*.RData | sed -e "s/RData/Rout/g")

.PHONY : analyses

analyses : $(simfiles)

%.Rout : %.RData
	# can pass args in by doing e.g. 'make ARGS="-l 3 -r 3" whatever.Rout'
	Rscript ising-inference.R -i $< -o $@ ${ARGS} 

test : check-ising.R
	Rscript --vanilla $<

info : 
	Rscript ../siminfo.R ising-sims
