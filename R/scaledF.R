scaled_F <- function(tml, eigenvalues){
  s1=sum(eigenvalues)
  s2=sum(eigenvalues^2)
  s3=sum(eigenvalues^3)
  
  denom<-2*s1*s2^2-s1^2*s3+2*s2*s3;
  if (denom>0) 
  {
    d1F3 <- s1*(s1^2*s2-2*s2^2+4*s1*s3)/denom;
    d2F3 <- (s1^2*s2 + 2*s2^2) / (s3*s1 - s2^2) +6; 
    if (d2F3<6) d2F3<-Inf;
    cF3 <- s1*(s1^2*s2 - 2*s2^2 + 4*s1*s3)/(s1^2*s2-4*s2^2+6*s1*s3); 
  }else
  {     
    d1F3<-Inf;
    d2F3<-s1^2/s2+4;
    cF3<- s1*(s1^2+2*s2)/(s1^2+4*s2);
  }
  unname(1-pf(tml/cF3,d1F3,d2F3))
}
