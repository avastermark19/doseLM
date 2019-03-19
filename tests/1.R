library(doseLM)
se <- simData(FA=100)
RUnit::checkEqualsNumeric(dim(se)[1]*dim(se)[2], 8000)
