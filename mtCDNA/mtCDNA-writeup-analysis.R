
library(contextutils)


load(gmfile)

mutpatinfo <- data.frame(
        mutpat = sapply( lapply( mutpats, lapply, paste, collapse="->" ), paste, collapse=" | " ),
        codon1 = factor( sapply( lapply( mutpats, lapply, "[", 1 ), paste, collapse=" | " ), levels=levels(codons$codon) ),
        codon2 = factor( sapply( lapply( mutpats, lapply, "[", 2 ), paste, collapse=" | " ), levels=levels(codons$codon) )
    )
mutpatinfo$aa1 <- codons$aa[ match(mutpatinfo$codon1,codons$codon) ]
mutpatinfo$aa2 <- codons$aa[ match(mutpatinfo$codon2,codons$codon) ]
mutpatinfo$synon <- with(mutpatinfo, aa1==aa2 )

for (k in seq_along(frame.mrun)) {
    mutpatinfo[ paste("frame",k,sep='') ] <- colMeans( frame.mrun[[k]]$batch[100:nrow(frame.mrun[[k]]$batch),1+1:length(mutpats)] )
}

with( mutpatinfo, cbind( levels(mutpat)[mutpat][order(frame1)], levels(mutpat)[mutpat][order(frame2)], levels(mutpat)[mutpat][order(frame3)] ) )
