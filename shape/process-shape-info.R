#!/usr/bin/Rscript

require(jsonlite)

# read in shape data from Jessica
x <- read.table("DNAshape-all-pentamers.txt",header=TRUE)

# rescale variables: chosen so SDs are nearish to 1.  Centering is not necessary, but might as well, to remind us the scale has changed.
xs <- scale(x[,-1], 
        center=c(
                HelT1=34.4,
                HelT2=34.4,
                MGW=5.0,
                ProT=-6.7,
                Roll1=-0.7,
                Roll2=-0.7
            ),
        scale=c(
                HelT1=1.5,
                HelT2=1.5,
                MGW=0.5,
                ProT=3.0,
                Roll1=3.0,
                Roll2=3.0
            ) )

# write out as JSON
y <- apply(xs,2,function (z) { names(z) <- x[,1]; as.list(z) })
cat(toJSON(y,pretty=TRUE),file="DNAshape-all-pentamers-scaled.json")
