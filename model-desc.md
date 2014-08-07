= Model description =

An object of class 'context' is the result of fitting a 'model' to some 'data'.
It should therefore carry around:
- 'data' : an object of class 'tuplecounts' giving numbers of observed paired tuples
- 'genmatrix' : an object of class 'genmatrix' with rows and columns indexed by 'headpats' (in the same order)
- 'mutrates' : instantaneous rates corresponding to 'genmatrix@mutpats'
- 'selcoef' : selection coefficients corresponding to 'genmatrix@selpats'
- 'params' : other parameters (branch lengths, arguments to fixfn, etc)
- 'projmatrix' : a (projection) matrix with rows indexed by 'headpats' and columns indexed by 'tailpats'
- 'likfun' : the function that returns the negative log-likelihood of the data as a function of (mutrates,selcoef,params)
- 'optim.results' : output from routine that estimated parameters (mutrates, selcoef, params)
Furthermore, it has the following methods:
- 'winlen( )' : an integer giving the length of the patterns that index rows of 'data'
- 'win( )' : an integer giving the length of the patterns that index columns of 'data'
- 'lwin( )' : an integer giving the offset that aligns columns of 'data' with rows of 'data'
- 'rownames( )' : a character vector of (winlen)-mers indexing the rows of 'data'
- 'colnames( )' : a character vector of (win)-mers indexing the columns of 'data'
- 'counts( )' : returns the (count) data
and even more usefully,
- 'coef( )' : a named list of mutation and selection coefficients
- 'fitted( )' : predicted counts under the model
- 'residuals( )' : fitted - counts
- 'residuals( , pretty=TRUE)' : residuals, in a data frame, with z-scores, sorted by z-score
Note that by default 'fitted(model)' predicts counts for the pattern lengths used to fit the model, but by passing other parameters, you can obtain predicted counts for other pattern lengths (but note: must pass in a new 'genmatrix').  The same goes for 'residuals(model)', but for patterns longer than those initially used, you must pass in the observed counts.

An object of class 'tuplecounts' is a matrix of counts of paired tuples, that additionally carries:
- 'counts' : a Matrix of counts
- 'lwin' : an integer giving the offset that aligns rownames(counts) with colnames(counts)


Note that although it is not stored in the 'context' object because it would be redundant, 'win' means the inner window width.

An object of class 'genmatrix' is a sparse matrix that additionally carries the following information (and some more stuff):
- 'bases' : character vector of allowed bases
- 'rownames( )' and 'colnames( )' : character vectors; should match with e.g. 'headpats' above.
- 'mutpats' : list of mutation motifs
- 'selpats' : list of selection motifs
- 'fixfn' : fixation function that translates differences in selection coefficient to mutation rate multipliers

