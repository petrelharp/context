#' @export
plot.contextMCMC <- function (model,...) {
    # plot the MCMC trace
    layout(t(1:2))
    batch <- model@results$batch
    params <- coef(model)
    matplot( batch[,1:nmuts(model)], type='l', lty=1:5, col=1:6 )
    abline(h=params[1:nmuts(model)], lty=4, lwd=2, col=1:6 )
    legend("topleft", legend=names(params)[1:nmuts(model)], lty=1:5, col=1:6, title="mutation" )
    matplot( batch[,nmuts(model)+(1:nsel(model))], type='l' )
    abline(h=params[nmuts(model)+(1:nsel(model))], lty=4, lwd=2, col=1:6 )
    legend("topleft", legend=names(params)[nmuts(model)+(1:nsel(model))], lty=1:5, col=1:6, title="selection" )
}

