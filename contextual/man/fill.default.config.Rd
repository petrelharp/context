% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input-output.R
\name{fill.default.config}
\alias{fill.default.config}
\title{Fill in Default Values in a Configuration List}
\usage{
fill.default.config(config, defaults = NULL)
}
\value{
A model configuration (list).
}
\description{
Fill in default values in a config list based on 'defaults', or length zeros if it makes sense
  for mutpats, mutrates, selpats, selfactors, selcoef, bases, fixfn, fixfn.params
}
\details{
Note that will NOT fill in default mutrates or selcoef if none are available in defaults
  and the corresponding patterns are not length zero
  (for instance, models for makegenmat don't need mutrates)
}
