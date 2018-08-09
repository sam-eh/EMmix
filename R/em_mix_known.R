#' EM Implementation for ancestry mixtures
#'
#' Implements the EM algorithm to estimate the mixture parameters of ancestries in an observed dataset using a reference dataset of some sort.
#'
#' @param x N by 3 data.frame for N markers and 3 genotype columns (e.g. "hom", "homref" and "het") giving allele frequencies of those observed at genetic marker.
#' @param MAF_thresh numeric vector of length 1 with value between 0 and 0.5 indicating when to exclude a marker due to rarity (currently not implemented).
#' @param pnames matrix of population names in column 1 that contains K rows (one for each ancestry considered), followed by the column names pertaining to the allele frequency(ies) for that population as they appear in dat.
#' @param pi_init numeric vector of length K providing initial mixture parameters.
#' @param threshold numeric vector of length 1 indicating when to terminate the loop. The loop will terminate when the sum of the absolute differences between the previous and current estimates for the population mixture is less than \code{threshold}.
#' @param dat data.frame object for frequencies by each population. This is the reference dataset that will be used.
#' @param Ntot total number of population considered for x.
#' @return EMmix object recording the iteration results from the function.
#' @examples
#' n <- 3000 #set number of markers to simulate
#' num_pop <- 5000 #set population size for observed data to be passed as x to em_mix_known
#' dat <- data.frame(AFR_AF = runif(n=n), NFE_AF = runif(n=n)) # simulate allele frequencies to be used for the reference dataset
#' l <- nrow(dat)
#' popNames <- c("AFR", "NFE")
#' hardyWeinNames <- paste(rep(popNames, each = 3), c("_hom","_het","_homref"),sep="")
#' dat_hardyWein <- as.data.frame(matrix(nrow=l, ncol=length(hardyWeinNames), dimnames = list(c(),hardyWeinNames)))
#' pnames <- matrix(nrow = length(popNames), ncol=4, dimnames = list(popNames, c()))
#' for(pop in popNames){
#'   column <- dat[,paste(pop, "_AF", sep = "")]
#'   # assume hardy-weinberg equilibrium
#'   hom <- column ^ 2
#'   het <- 2 * column * (1-column)
#'   homref <- (1-column)^2
#'   cols <- paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")
#'   dat_hardyWein[,cols] <- cbind(hom, het, homref)
#'   pnames[pop,] <- c(pop, cols)
#' }
#' real_pi <- c(0.4,0.6) #set the real proportions expected for the simulated observed data
#' pop_number <- n  * real_pi
#' names(pop_number) <- popNames
#' names(real_pi) <- popNames
#' pop_sim<-as.data.frame(matrix(nrow = l, ncol = 0))
#' for(pop in popNames){
#'   popDat <- dat_hardyWein[,paste(rep(pop, each = 3), c("_hom","_het","_homref"),sep="")]
#'   pop_sim <<- cbind(pop_sim, t(apply(popDat, 1, function(pd){rmultinom(n=1, size = pop_number[pop], prob = pd)})))
#' }
#' # next need to get the overall frequency from the simulated columns
#' sim_x <- data.frame(hom = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 1)], 1, sum),
#'                     het = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 2)], 1, sum),
#'                     homref = apply(pop_sim[,which(1:(3 * length(popNames)) %% 3 == 0)], 1, sum)
#' )
#' mixture <- EMmix::em_mix_known(x = sim_x, dat = dat_hardyWein, pnames = pnames, pi_init = c(0.1,0.9), Ntot = NA)
#' #' @export
em_mix_known <- function(x, dat, pnames, pi_init, Ntot, threshold = 0.01, MAF_thresh = 0.05, path = NA){
  pi_out <- pi_out_median <- pi_out_90 <- pi_new <- pi_median <- pi_90 <- pi <- pi_init
  names(pi_median) <- names(pi_new) <- names(pi_90) <- names(pi) <- pnames[,1]
  k <- nrow(pnames)
  N <- nrow(x)
  thresh_check <- threshold + 1
  cols_toSelect <- character(); apply(as.matrix(pnames[,2:ncol(pnames)]), 1, function(pop){cols_toSelect <<- c(cols_toSelect, pop)} )
  p <- dat # p <- dat(which(x$AF > MAF_thresh & x$AF < (1-MAF_thresh)))
  iter = 0
  if(ncol(pnames) < 2){stop("something is wrong with the dimensions of the value provided for pnames")}
  if(ncol(pnames)>2){
    while(thresh_check > threshold){
      pi = pi_new
      # expectation step
      gamma_tmp.a <- cbind(
        apply(pnames, 1, function(pop){
          x.p <- cbind(x, p[,pop[2:ncol(pnames)]])
          apply(x.p, 1, function(gen){
            pi[pop[1]]*dmultinom(as.numeric(gen[1:(ncol(pnames)-1)]), prob = gen[ncol(pnames):(2 * (ncol(pnames)-1))], log = T)
          })
        })
      )
      gamma_tmp<-gamma_tmp.a/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})

      # maximization step
      pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
      names(pi_new) <- names(pi)
      pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
      names(pi_median) <- names(pi)
      pi_90<-apply(gamma_tmp,2,function(x){quantile(x, probs=0.90,na.rm=T)})
      names(pi_90) <- names(pi)

      pi_out<-rbind(pi_out, pi_new)
      pi_out_median<-rbind(pi_out_median, pi_median)
      pi_out_90<-rbind(pi_out_90, pi_90)
      iter=iter+1
      thresh_check<-sum(abs(pi-pi_new))
    }
  }
  # This accounts for instances when x is N x 1 (and pnames is therefore only K x 2); x is N x 1 because the one column is allele frequency
  if(ncol(pnames) == 2){
    if(any(x[,1] < 1)){x[,1] <- round(x[,1] * 2 * Ntot)} # need x data as counts for dbinom function
    while(thresh_check > threshold){
      pi = pi_new
      # expectation step
      gamma_tmp.a <- cbind(
        apply(pnames, 1, function(pop){
          x.p <- cbind(x, p[,pop[2:ncol(pnames)]])
          pi[pop[1]]*dbinom(x = x.p[,1], size = Ntot, prob = x.p[,2], log = TRUE) # this only works when x is given as integers and not frequencies
        })
      )
      gamma_tmp<-gamma_tmp.a/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})

      # maximization step
      pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
      names(pi_new) <- names(pi)
      pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
      pi_90<-apply(gamma_tmp,2,function(x){quantile(x, probs=0.90,na.rm=T)})

      pi_out<-rbind(pi_out, pi_new)
      pi_out_median<-rbind(pi_out_median, pi_median)
      pi_out_90<-rbind(pi_out_90, pi_90)
      iter=iter+1
      thresh_check<-sum(abs(pi-pi_new))
    }
  }
  colnames(gamma_tmp) <- pnames[,1]
  mix <- EMmix(pnames = pnames, x = as.data.frame(x), dat = as.data.frame(dat),
               pi_init = pi_init, med = pi_out_median, mu = pi_out, x90 = pi_out_90, geneMixes = gamma_tmp)
  show(mix)
  return(mix)
}
