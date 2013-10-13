#/usr/bin/R --vanilla

source("../sim-context-fns.R")

bases <- c("X","O")
mutpats <- list( 
    list( c("X","O"), c("O","X"), c("XO","OX") )
    ) 
mutrates <- 1
selpats <- list()
selcoef <- numeric(0)
