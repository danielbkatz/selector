
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


eigen_block <- function(eigenvalues, numblocks){
  if (numblocks <2)
    return(NA)
  e = eigenvalues
  t=chunk(eigenvalues, numblocks)
  m = sapply(t, mean)
  counter =1
  for(i in 1:length(t)){
    len = length(t[[i]])
    e[counter:(counter+len-1)] <- m[i]
    counter <- counter + len
  }
  return(e)
}
