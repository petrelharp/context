simfiles = $(wildcard ising-sims/*.RData)
outfiles = $(patsubst %.RData,%.Rout,$(simfiles))

.PHONY : analyses mcmcs

analyses : $(outfiles)

%.Rout : %.RData
	# can pass args in by doing e.g. 'make ARGS="-l 3 -r 3" whatever.Rout'
	Rscript ising-inference.R -i $< -o $@ ${ARGS} 

mcmcs : 
	# run one more batch of mcmcs for everything
	for x in $(simfiles); do Rscript ising-mcmc.R -i $$x -n 10000; done

test : check-ising.R
	Rscript --vanilla $<

info : 
	Rscript ../siminfo.R ising-sims
