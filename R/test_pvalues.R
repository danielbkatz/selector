#should be public
test_pvalues <- function(f){
  
  UG <- lavInspect(f, "UGamma")
  df <- fitmeasures(f, "df")
  tml <- fitmeasures(f, "chisq")
  
  porigs <- c(ml = unname(fitmeasures(f, "pvalue")), 
              get_pvalues(tml,  Re(eigen(UG)$values[1:df])),
              ss = unname(get_ss_pvalue(tml, df, UG)))
  porigs
}