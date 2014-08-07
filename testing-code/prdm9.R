

bases <- c("A","T","C","G")

source("../sim-context-fns.R")
source("../context-inference-fns.R")

# maximum size of pattern (for simulation)
mutpats <- c(
        apply(combn(bases,2),2,list),  # single-base rates
        apply(combn(bases,2)[2:1,],2,list),  # single-base rates
        list( list( c("CG","TG"), c("CG","CA") ) )  # CpG rate
    ) 
mutrates <- c( rep(2e-8,length(mutpats)-1), 2e-7 )
selpats <- list("CCGCCGT")
selcoef <- 1e-4
fixfn <- popgen.fixfn

initfreqs <- rep( 1/4, 4 )

seqlen <- 172
initseq <- rinitseq(seqlen,bases,basefreqs=initfreqs)

system.time( 
        simseqs <- list(
                simseq( seqlen=seqlen, tlen=tlen, mutpats=mutpats, mutrates=mutrates, selpats=selpats, selcoef=selcoef, initseq=initseq, bases=bases )
            )
    )

###
# testing

selpats <- list("CCGCC")
longwin <- 8
genmatrix <- makegenmatrix( patlen=longwin, mutpats=mutpats, selpats=selpats, boundary='none', Ne=1e4 )
