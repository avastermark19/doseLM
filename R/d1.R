#' @title d1 Function to run Tom's ANOVA on se object.
#' @description This function runs Tom's ANOVA on se object.
#' @usage d1(se)
#' @param se SummExp object.
#' @details This function runs Tom's ANOVA on se object.
#' @return Returns/prints Tom's ANOVA.
#' @examples
#' type1 <- simData()
#' d1(type1)
#' @author AJ Vaestermark, JR Walters.
#' @references The 'doseLM' package, 2019 (in press).

d1 <- function(se) {

#library(doseR)

  #data(hmel.se)
  f_se <- quantFilter(se, lo.bound = 0.35, hi.bound = 0.65)
  dm <- se.DM(f_se)
  #glSeq(dm, "-1 + replicate")
    #dm<- cD.DM(cd.leg)
  gl  <-  tryCatch(glSeq(dm, "-1 + replicate"), error=function(e) NA)
  #gl2 <- tryCatch(glSeq(dm, "-1 + annotation.ZA + replicate"), error=function(e) NA)
  #gl3 <- tryCatch(glSeq(dm, "-1 + annotation.ZA * replicate"), error=function(e) NA)

  gl2 <- tryCatch(glSeq(dm, "-1 + annotation.something + replicate"), error=function(e) NA)
  gl3 <- tryCatch(glSeq(dm, "-1 + annotation.something * replicate"), error=function(e) NA)
  gl4 <- tryCatch(glSeq(dm, "-1 + annotation.something"), error=function(e) NA)

  # "-1 + replicate*annotation.ZA"

  print(   anova(gl, gl2)[8]   )
  #print ( "------------------------------------------------")
  print(   anova(gl, gl3)[8]   )
  #print ( "------------------------------------------------")
  print(   anova(gl, gl4)[8]   )
  #print ( "------------------------------------------------")

  print(   anova(gl2, gl3)[8]   )
  #print ( "------------------------------------------------")
  print(   anova(gl2, gl4)[8]   )
  #print ( "------------------------------------------------")

  print(   anova(gl3, gl4)[8]   )
  #print ( "------------------------------------------------")



}


# Type 1
# type1 <- simData() ; # you get sig result for the comparison between just MF and just AZ
# Type 2
# type2 <- simData(FZ=1, MZ=1) ; you get sig result for comparison between just MF and all the other 3
# Type 3
# type3 <- simData(MZ=1) ; # you get sig result for all comparisons except: justMF and MF+AZ, just MF and just AZ
# Type 4
# type4 <- simData(FZ=100) ; # interaction effect, secondd one sig (same as type 3)
# "TYpe 5"
# type5 <- simData(FA= 100, FZ=100) ; you get sig result for all comparisons except: justMF and MF+AZ, just MF and just AZ


#       MF	  MF	  MF	MF+AZ	MF+AZ	MF*AZ
#       MF+AZ	MF*AZ	AZ	MF*AZ	AZ	  AZ
#Type1  1	    1	    ***	1	    1	    1
#Type2  ***	  **	  ***	1	    1	    1
#Type3  1	    ***	  1	  ***	  ***	  ***
#Type4  *	    ***	  1	  ***	  ***	  ***
#Type5  1	    *	    1	  *	    ***	  ***

# se <- simData(FA=100) ; 0.47 and 2.2e-16
