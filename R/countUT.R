#' npc
#'
#' Function to calculate the npc
#'
#' See Ma et al. (2019) for all the details.
#'
#'
#' @param output Must be an array storing the probabilities attributed to \code{nobs}
#' patients to the 3 different benefit levels of the two competing treatments
#' for \code{nset} replicas of the simulation study.
#' Its dimension are \code{nobs} x 6 x \code{nset}, where the first three columns are
#' the probabilities for treatment 1 and the latter
#' 3 columns store the probabilities for treatment 2.
#' @param trtsgn treatment assigned to the patients by design.
#' @param myoutot List of \code{nset} outcome variables in ordinal notation
#'
#' @references
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' @return a \code{nobs} vector storing the npc for each patient,
#' @export

npc <- function(output, trtsgn, myoutot){
   K <- dim(output)[3]
   n <- dim(output)[1]
   myctut <- array(0, dim = c(3, 3, K))
   myctutSum <- NULL

   for(i in 1:K){
      mycurdata <- output[,,i]
      mypre <- NULL
      pretrt1 <- apply(mycurdata[,1:3], 1, which.max)
      pretrt2 <- apply(mycurdata[,4:6], 1, which.max)
      mypreTall <- cbind(pretrt1, pretrt2)

      for(j in 1:n){
         mypre[j] <- mypreTall[j, trtsgn[j]]
      }
      sts <- table(mypre, myoutot[,i])
      mysdls <- as.numeric(rownames(sts))
      str1 <- matrix(0, nrow = 3, ncol = 3)
      str1[mysdls,] <- sts

      myctut[,,i] <- str1*diag(3)
      myctutSum[i] <- sum(str1*diag(3))
   }

   res <- cbind(myctutSum)
   return(res)
}
