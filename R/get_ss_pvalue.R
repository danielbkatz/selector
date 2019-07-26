get_ss_pvalue <- function(tml, d, UG){
  
  tr_UG2 <- sum(sapply(X=1:ncol(UG), FUN=function(i) sum(UG[i,]*UG[,i]))) 
  trUG <- tr(UG)
  c1 <- sqrt(d/tr_UG2)
  c2 <- d-sqrt(d*trUG^2/tr_UG2)
  
  1-pchisq(c1*tml+c2, df=d)
}