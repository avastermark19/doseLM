#' @title d0 Function to run limma voom edgeR on se object.
#' @description This function runs limma voom edgeR on se object.
#' @usage d0(se)
#' @param se SummExp object.
#' @details This function runs limma voom edgeR on se object.
#' @return Returns/prints voom plots.
#' @examples
#' type1 <- simData()
#' d0(type1)
#' @author AJ Vaestermark, JR Walters.
#' @references The 'doseLM' package, 2019 (in press).

d0 <- function(se) {

#library(edgeR)

#library(doseR)

 # data(hmel.data.doser)

#counts <- as.data.frame(hmel.dat$readcounts)
  counts <- as.data.frame(   assays(se)$counts )

d0 <- DGEList(counts)

d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left

#group <- as.factor(c("Male", "Male", "Male", "Female", "Female", "Female"))
group <- as.factor(c("Female", "Female", "Male", "Male"))

limma::plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0 + group)

y <- suppressWarnings(limma::voom(d, mm, plot = T))

tmp <- suppressWarnings(limma::voom(d0, mm, plot = T))

fit <- limma::lmFit(y, mm)
head(coef(fit))

groupMale <- NULL
groupFemale <- NULL

contr <- limma::makeContrasts(groupMale - groupFemale, levels = colnames(coef(fit)))

tmp <- limma::contrasts.fit(fit, contr)

tmp <- limma::eBayes(tmp)

top.table <- limma::topTable(tmp, sort.by = "P", n = Inf)

print(head(top.table, 20))

#length(which(top.table$adj.P.Val < 0.05))

}
