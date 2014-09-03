library("pander")

load("sim-tasep-123456-genmatrix-4-complete-54321.RData")

# le sigh.
n_hashes <- function(n) if(n > 0) paste0("#", n_hashes(n-1))

custom_pander <- function(x, sec_depth=1) {
    x_class = class(x)
    # first catch any S4 objects that want special treatment
    if (x_class == "dgTMatrix" || x_class == "dgCMatrix") {
        pander(as.matrix(x))
    }
    # S4 backup
    else if (isS4(x)) {
        sapply(
            names(attributes(x)),
            function(name) {
                y <- slot(x, name)
                cat(n_hashes(sec_depth), name, ":", class(y), "\n")
                custom_pander(y, sec_depth+1)
            })
    }
    else if(x_class == "function") { }
    else if (x_class == "list") {
        if(length(x) > 0) pander(x)
    }
    else {
        if(length(x) > 0) cat(x)
    }
    cat("\n\n")
}
s <- capture.output(custom_pander(model))

cat(s,sep="\n")

