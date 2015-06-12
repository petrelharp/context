#!/usr/bin/Rscript

source("../context-inference-fns.R",chdir=TRUE)
source("../sim-context-fns.R",chdir=TRUE)
source("../input-output.R",chdir=TRUE)

# this has X's as VERY deleterious; O's as VERY beneficial
#  and  sp2 on a LONG branch
config <- parse.models( treeify.config( read.config("varying-selcoef.json") ) )

simseqs <- simseq.tree(100,config)

longpats <- getpatterns(4,bases=config$bases)
shortpats <- getpatterns(1,bases=config$bases)

# counts
shorts.1 <- counttrans(shortpats, shortpats, simseqs=simseqs$sp1 )
shorts.2 <- counttrans(shortpats, shortpats, simseqs=simseqs$sp2 )

# for sp1, should have no O->X
stopifnot( ( counts(shorts.1)["O","X"] == 0 ) )
# for sp2, should have no X's at all
stopifnot( ( counts(shorts.2)["O","X"] == 0 ) && ( counts(shorts.2)["X","X"] == 0 ) )


