# Iterative procedure for producing variance estimates
test.cov <- function(n1,n2,m1,m2,s1,s2){
  # n = sample size [scalar]
  # m = sample mean in ilr coordinates [1,D-1]
  # s = covariance matrices of two sample populations
  source("common_mean.R")
  # Initialize ratios != 1
  ratio1 = 0
  ratio2 = 0
  # # Initialize a counter to keep track of the number of iterations needed for convergence
  h = 0
  # # Initialize intermediate covariance matrices using the initial values; they will be reset during iteration
  sa = s1
  sb = s2
  while (TRUE){
    h = h+1
    mhi = common.mean(n1,n2,sa,sb,m1,m2)
    sc1 = s1+(t(m1-mhi)%*%(m1-mhi))
    sc2 = s2+(t(m2-mhi)%*%(m2-mhi))
    ratio1 = round(det(sa)/det(sc1),6)
    ratio2 = round(det(sb)/det(sc2),6)
    if (ratio1 && ratio2 == 1){
      test_c = list(sc1,sc2)
      break
    } else {
      sa = sc1
      sb = sc2
    }
  }
  print(h)
  return(test_c)
}