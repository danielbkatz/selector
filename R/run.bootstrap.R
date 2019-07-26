#worker function
# X is seed, startvals is a parTable object
run_bootstrap <- function(X, b.reps, sample.transformed, startvalues){
  #possible statistics
  set.seed(X)
  pb <- txtProgressBar(min = 0, max = b.reps, style = 3)
  result.df <- data.frame(matrix(nrow=b.reps, ncol=length(porigs)))
  for(b in 1:b.reps){
    setTxtProgressBar(pb, b)
    boot.sample <- sample.transformed[sample(1:nrow(sample.transformed), replace=T), ]
    
    f.boot <- tryCatch(sem(startvalues,boot.sample, start=startvalues), error=function(w) { TRUE})
    while(isTRUE(f.boot)){
      boot.sample <- sample.transformed[sample(1:nrow(sample.transformed), replace=T), ]
      f.boot <- tryCatch(sem(startvalues,boot.sample, start=startvalues), error=function(w) { TRUE})
    }
    
    UG.boot <- lavInspect(f.boot, "UGamma")
    tml.boot <- fitmeasures(f.boot, "chisq")
    
   result.df[b, ] <- c(ml = unname(fitmeasures(f.boot, "pvalue")), 
                get_pvalues(tml.boot,  Re(eigen(UG.boot)$values[1:df])),
                ss = unname(get_ss_pvalue(tml.boot, df, UG.boot)))
  }
  return(result.df)
}

