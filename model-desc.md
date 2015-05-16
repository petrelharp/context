= Model description =

An object of class 'context' is the result of fitting a 'model' to some 'counts' data.
It should therefore carry around:
- 'counts' : an object of class 'tuplecounts' giving numbers of observed paired tuples
- 'genmatrix' : an object of class 'genmatrix' with rows and columns indexed by 'headpats' (in the same order)
- 'mutrates' : instantaneous rates corresponding to 'genmatrix@mutpats'
- 'selcoef' : selection coefficients corresponding to 'genmatrix@selpats'
- 'params' : other parameters (branch lengths, arguments to fixfn, etc)
- 'projmatrix' : a (projection) matrix with rows indexed by 'headpats' and columns indexed by 'tailpats'
- 'likfun' : the function that returns the negative log-likelihood of the data as a function of (mutrates,selcoef,params)
- 'results' : output from routine that estimated parameters (either optim() or metrop())
- 'invocation' : an optional character string giving the command that produced it
Furthermore, it has the following methods:
- 'winlen( )' : an integer giving the length of the patterns that index rows of 'counts'
- 'win( )' : an integer giving the length of the patterns that index columns of 'counts'
- 'lwin( )' : an integer giving the offset that aligns columns of 'counts' with rows of 'counts'
- 'rownames( )' : a character vector of (winlen)-mers indexing the rows of 'counts'
- 'colnames( )' : a character vector of (win)-mers indexing the columns of 'counts'
and even more usefully,
- 'coef( )' : a named list of mutation and selection coefficients (mutrates, selcoef, params)
- 'fitted( )' : predicted counts under the model
- 'residuals( )' : fitted - counts
- 'residuals( , pretty=TRUE)' : residuals, in a data frame, with z-scores, sorted by z-score
Note that by default 'fitted(model)' predicts counts for the pattern lengths used to fit the model, but by passing other parameters, you can obtain predicted counts for other pattern lengths (but note: must pass in a new 'genmatrix').  The same goes for 'residuals(model)', but for patterns longer than those initially used, you must pass in the observed counts.

An object of class 'tuplecounts' is a matrix of counts of paired tuples, that additionally carries:
- 'lwin' : an integer giving the offset that aligns long patterns with short patterns
- 'bases' : a character vector of allowed bases
- 'counts' : a Matrix of counts, with rows indexed by long patterns and columns indexed by (combinations of) short patterns; and
  - 'rownames(counts)' : long patterns.
  - 'colnames(counts)' : either short patterns (for two taxa), or arbitrary, with corresponding information encoded by colpatterns
- 'rowtaxon' : the taxon in which we count long patterns
- 'colpatterns' : a data frame whose columns are indexed by taxa in which short patterns are counted, and the k-th row of which gives the combination of short patterns corresponding to the k-th column of counts
- 'coltaxa( )' : the column names of colpatterns.
For instance, if the tree has taxa `sp1`, `sp2`, and `sp3`; `bases` is `A,T`; and `rowtaxon` is `sp1`; then `rownames(counts)` could be `AA,AT,TA,TT`, and the rows of `colpatterns` could be `A,A`, `A,T`, `T,A`, and `T,T`; meaning that, for instance, there were `counts[1,2]` times that `sp1` was found to have `AA` at a site where `sp2` had `A` while `sp3` had `T`.
Methods:
- 'countframe' : returns the `counts` as a data.frame, with the first, named, columns giving the patterns in each taxa, and the last (named `count`) giving the number of occurrences


An object of class 'genmatrix' is a sparse matrix that additionally carries the following information (and some more stuff):
- 'bases' : character vector of allowed bases
- 'rownames( )' and 'colnames( )' : character vectors; should match with e.g. 'headpats' above.
- 'mutpats' : list of mutation motifs
- 'selpats' : list of selection motifs
- 'fixfn' : fixation function that translates differences in selection coefficient to mutation rate multipliers
- 'nmuts( )': the number of mutation patterns present in that genmatrix


= Config files =

A full config file can have:
- 'tree' : tree, in Newick format, with node labels and optional edge lengths, 
- 'bases' : as above
- 'named model stanzas' : including one for each node label, specifying the model that occurs on the branch above the named node.  This can also be a character string referring to a different named model stanza, indicating that the two edges should have the same model, and 'share parameters'.
Other stuff (e.g. "comment") will be ignored.

Tree-based config files should also have:
- 'initfreqs' : base frequencies at the root.
- 'initfreqs.scale' : scaling factor for initfreqs.
- 'tlen.scale' : scaling factor for branch lengths in the tree.

A 'model stanza' has:
- 'mutpats' : list of lists of character pairs (one list per mutation motif)
- 'mutrates' : numeric, nonnegative, same length as mutpats
- 'selpats' : either:
  - list of lists of character strings (one list per selection motif)
    - *and* [optionally] : `selfactors' : list of numeric vectors
  - **or** : list of lists of named numeric vectors, one list per selection motif, and the numbers multiplying the selection coefficient (names will be used as selpats; numbers as selfactors)
- 'selcoef' : numeric, same length as selpats
- 'fixfn' : name of a function, or R code (e.g. "function (x) { ... }"
- 'fixfn.params' : named list of additional parameters to fixfn
- 'genmatrix' : pattern for where to save genmatrix files, with `%` to be substituted for the pattern length
- 'mutrates.prior', 'selcoef.prior', 'fixfn.params.prior' : coefficients for priors on respective parameters
- 'mutrates.scale', 'selcoef.scale', 'fixfn.params.scale' : scale over which to move in optimization, MCMC, etcetera. 

Parameters whose scale is set to zero *will be held as fixed* in the analysis.
