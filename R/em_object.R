#' EMmix Object
#'
#' Records information from the EM iterations when running em_mix_known
#'
#' @param pnames a matrix of ancestry names in the first column and the corresponding column names for each genotype.
#' @param x the observed data.
#' @param dat the reference data that is used to determine the mixture proportions in x.
#' @param pi_init the initiating mixture proportions.
#' @param med the median values of the estimated mixture proportions of each ancestry from each EM iteration.
#' @param mu the mean values of the estimated mixture proportions of each ancestry from each EM iteration.
#' @param x90 the upper 90th quantile values of the estimated mixture proportions of each ancestry from each EM iteration.
#' @param geneMixes the final estimates for the genetic inheritances. Each row indicates a different gene, and each column indicates an ancestry.
#' @export

EMmix <- setClass(
  "EMmix",
  slots = c(
    pnames = "matrix", x="data.frame", dat="data.frame", pi_init="numeric",
    med = "matrix", mu = "matrix", x90 = "matrix", geneMixes = "matrix"
  )
)

setMethod("initialize",
          signature = c(.Object="EMmix"),
          function(.Object, pnames, x, dat, pi_init, med, mu, x90, geneMixes){
            .Object@pnames = pnames
            .Object@x = x
            .Object@dat = dat
            .Object@pi_init = pi_init
            .Object@med = med
            .Object@mu = mu
            .Object@x90 = x90
            .Object@geneMixes = geneMixes
            return(.Object)
          }
)

#' Plot EM Object
#'
#' Plots the iteration results to display the rate of convergence of the algorithm to the final results
#'
#' @param x an EMmix object
#' @export

setMethod("plot",
          signature = c(x="EMmix",y="missing"),
          function(x,y){
            last_iter <- nrow(x@med)
            pops <- x@pnames[,1]
            for(p in pops){
              xmin <- min(x@med[,p])
              xmax <- max(x@med[,p])
              last <- x@med[last_iter,p]
              plot.default(x=x@med[,p], main=p, xlab="Iteration Number", ylab = "Median Mixture Parameter", ylim=c(xmin,xmax),type="p")
              abline(h=last, lty=2, col="firebrick1")
            }
          }
)

#' Show EMmix Object
#'
#' Show function for an EM mix object
#'
#' @param object an EMmix object
#' @export

setMethod("show",
          signature = "EMmix",
          function(object){
            last_iter <- nrow(object@med)
            cat(last_iter," iterations conducted\n\nFinal Median Mixture Parameters:\n")
            print(sort(object@med[last_iter,], decreasing = TRUE))
            cat("\n\nFinal Mean Mixture Parameters:\n")
            print(sort(object@mu[last_iter,], decreasing = TRUE))
          })
