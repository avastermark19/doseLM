#' @title simData Function to make simData.
#' @description This function generates simData.
#' @usage simData(FA, MA, FZ, MZ)
#' @param FA FA fold change parameter.
#' @param MA MA fold change parameter.
#' @param FZ FZ fold change parameter.
#' @param MZ MZ fold change parameter.
#' @details This function makes simData.
#' @return Returns simData.
#' @examples
#' se <- simData(FA=100)
#' str(se)
#' @author AJ Vaestermark, JR Walters.
#' @references The 'doseLM' package, 2019 (in press).

simData <- function(FA = 10, MA = 10, FZ = 10, MZ = 10) {

  samples <- c('F','F','M','M')
  something <- rep("A", 2000)
  something[seq(1000,2000)]  <- "Z" ###########################################################################

  FC <- matrix(nrow = length( levels(factor(something))  )  , ncol = length(levels(factor(samples))))

  colnames(FC) <- levels(factor(samples))
  rownames(FC) <- levels(factor(something))

  FC["A","F"] <- FA
  FC["A","M"] <- MA
  FC["Z","F"] <- FZ
  FC["Z","M"] <- MZ

  Ngenes <- 2000
  minLength <- 55
  trxL <- rgamma(Ngenes, shape = 1.5, scale = 1000)
  trxL[trxL < minLength] <- minLength

  trxN <- rgamma(Ngenes, shape = .5, scale = 1)
  readspertx <- round(trxL * trxN) + 1

  size.scale <- 5

  Nreps <- 4

  ctsmtx <- matrix(nrow = Ngenes, ncol = length(samples) )
  dim(ctsmtx)

  for (i in seq_len(length(samples))){ #################################################################
    for (j in seq_len(Ngenes)){ #######################################################################
      ctsmtx[j,i] <- rnbinom(n = 1, mu = readspertx[j]  * FC[ something[j] ,
      samples[i] ], size = readspertx[j]/size.scale )
    }
  }

  ###########################

  counts <- ctsmtx
  colData <- S4Vectors::DataFrame(Treatment=as.factor(samples) ,
                                  row.names=seq_len(4)) #######################################################
  rowData <- S4Vectors::DataFrame(annotation = as.data.frame(something),
                                  seglens = trxL,
                                  row.names=seq_len(2000) ) ###########################################################
  se3 <- SummarizedExperiment(assays=list(counts=counts),
                              colData=colData, rowData=rowData)

  SummarizedExperiment::colData(se3)$Libsizes <- c(2000,2000,2000,2000)
  SummarizedExperiment::assays(se3)$rpkm <- make_RPKM(se3)

  plotExpr(se3, groupings = "annotation.something", clusterby_grouping = FALSE,
           col=c("red","red","blue","blue"), notch=TRUE, ylab = "Log2(RPKM)", xlab='')
  invisible(se3)
}
