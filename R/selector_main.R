rm(list=ls())
#
#install these packages first, if necessary
#
library(psych)
library(CompQuadForm)
library(reshape2)
library(lavaan)
library(ggplot2)
library(goftest)

#worker function. Easy to put in parallel
run.bootstrap <- function(X, b.reps, sample.transformed, startvals){
  set.seed(X)
  pb <- txtProgressBar(min = 0, max = b.reps, style = 3)
  pml=NULL; psb=NULL; pss=NULL;
  for(b in 1:b.reps){
    setTxtProgressBar(pb, b)
    boot.sample = sample.transformed[sample(1:nrow(sample.transformed), replace=T), ]
    f.boot = tryCatch(sem(model,boot.sample, start=startvals), error=function(w) { TRUE})
    f.boot.sb = tryCatch(sem(model, boot.sample, test="satorra.bentler", start=parTable(f.boot)),error=function(w) { TRUE})
    f.boot.ss = tryCatch(sem(model, boot.sample, test="scaled.shifted", start=parTable(f.boot)),error=function(w) { TRUE})

    pml <- c(pml,fitmeasures(f.boot,"pvalue"))
    psb <- c(psb,fitmeasures(f.boot.sb,"pvalue.scaled"))
    pss <- c(pss,fitmeasures(f.boot.ss,"pvalue.scaled"))

  }
  return(data.frame(pml,psb,pss))
}



selector <- function(f, b.reps=1000, seed=1){
  #continuous case
  f.orig <- f
  orig.sample <- lavInspect(f.orig, what="data")
  model <- parTable(f.orig)

  porig <- fitmeasures(sem(model, orig.sample), "pvalue")
  porig <- c(porig, fitmeasures(sem(model, orig.sample, test="satorra.bentler"), "pvalue.scaled"))
  porig <- c(porig, fitmeasures(sem(model, orig.sample, test="scaled.shifted"), "pvalue.scaled"))


  #####
  # transform sample a la Bollen-Stine
  Sigma.hat <- lavaan:::computeSigmaHat(lavmodel = f.orig@Model)
  sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[1]])
  S.inv.sqrt <- lav_matrix_symmetric_sqrt(f.orig@SampleStats@icov[[1]])
  sample.transformed <- data.frame(as.matrix(orig.sample) %*% S.inv.sqrt %*% sigma.sqrt)
  colnames(sample.transformed) <- colnames(orig.sample)
  #####

  ## run bootstrap in serial
  res.df <- run.bootstrap(X=seed, b.reps=b.reps, sample.transformed=sample.transformed, startvals=model)
  colnames(res.df)<- c("ML", "SB", "SS")

  distUniform= as.vector(sapply(res.df,function(x) ad.test(x, null="punif")$statistic ))

  selected <- which.min(distUniform)
  selectedname <- colnames(res.df)[selected]
  pselected <- porig[selected]

  cat("The selector selects: ", selectedname, " with p-value", pselected)

  rownames(res.df) = NULL
  melted = melt(res.df)
  melted$statistic <- melted$variable

  #plot of the p-values
  ggplot(melted, aes(value))+geom_histogram()+facet_wrap(~statistic)

}



