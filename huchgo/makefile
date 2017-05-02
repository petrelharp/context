simfiles = $(wildcard huchgo-sims/*.RData)
outfiles = $(patsubst %.RData,%.Rout,$(simfiles))
# mcmcfiles = $(shell for $$x in $(patsubst %.Rdata,,$(simfiles)); do \
# 			k=0; for ii in $(find $$x -name *mcmc*.Rout); do (( k++ )); done\
# 			echo $$x/mcmc-run-$$k.Rout; done)
# mcmcfiles = $(shell for x in $(patsubst %.RData,%,$(simfiles)); do k=0; for ii in $(find $$x -name "*mcmc*.Rout"); do let "k+=1"; done; echo $$x/mcmc-run-$$k.Rout; done)
mcmcfiles = $(shell for x in $(patsubst %.RData,%,$(simfiles)); do k=$$(if [ -e $$x ]; then ls $$x | grep mcmc-run | grep Rout | wc -l; else echo 0; fi); echo $$x/mcmc-run-$$((k+1)).Rout; done)

.PHONY : analyses mcmcs

analyses : $(outfiles)

mcmcs : $(mcmcfiles)

%.Rout : %.RData
	# can pass args in by doing e.g. 'make ARGS="-l 3 -r 3" whatever.Rout'
	Rscript huchgo-inference.R -i $< -o $@ ${ARGS} 

info : 
	Rscript siminfo.R huchgo-sims > $@
	Rscript parse-mcmc.R  >> $@
	Rscript pairplots-mcmc.R huchgo-sims

%/mcmc-run-1.Rout : %.RData
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-2.Rout : %.RData %/mcmc-run-1.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-3.Rout : %.RData %/mcmc-run-2.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-4.Rout : %.RData %/mcmc-run-3.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-5.Rout : %.RData %/mcmc-run-4.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-6.Rout : %.RData %/mcmc-run-5.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-7.Rout : %.RData %/mcmc-run-6.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-8.Rout : %.RData %/mcmc-run-7.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-9.Rout : %.RData %/mcmc-run-8.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}

%/mcmc-run-10.Rout : %.RData %/mcmc-run-9.Rout 
	Rscript huchgo-mcmc.R -i $< -o $@ ${ARGS}
