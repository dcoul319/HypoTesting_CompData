# Hypothesis Testing Function
# Follows the procedure outlined in Chapter 7.3 of 'Modeling and Analysis of Compositional Data' (2015)
# by Pawlowsky-Glahn, Egozcue, and Tolosana-Delgado
ml.ratiotest <- function(n1,n2,m1,m2,s1,s2){
  # n = sample size [scalar]
  # m = sample mean in ilr coordinates [1,D-1]
  # s = covariance matrices of two sample populations
  D = length(m1)+1 # Used to calculate the degrees of freedom in the ChiSquare Distribution
  # Call functions for the iterative procedure used in the third hypothesis test
  source("common_mean.R")
  source("test_covariance.R")
  # Calculate the pooled covariance matrix as well as the combined sample estimates
  pool_cov = ((n1*s1)+(n2*s2))/(n1+n2)
  comb_mean = ((n1*m1)+(n2*m2))/(n1+n2)
  comb_cov = pool_cov+(n1*n2*(m1-m2)%*%t(m1-m2))/((n1+n2)^2)
  # Calculate the test statistics for the first two hypothesis tests
  # First test: same center and covariance structure vs general hypothesis
  Qavg = n1*log(det(comb_cov)/det(s1))+n2*log(det(comb_cov)/det(s2))
  deg_free = 0.5*(D-1)*(D+2)
  CHIavg = qchisq(0.67,deg_free,lower.tail=FALSE)
  if (Qavg >= CHIavg){
    print("H1:FAILURE")
  } else {
    print("H1:PASS")
  }
  print(deg_free)
  # Second test: different center but same covariance structure
  Qbvg = n1*log(det(pool_cov)/det(s1))+n2*log(det(pool_cov)/det(s2))
  deg_free = 0.5*D*(D-1)
  CHIbvg = qchisq(0.67,deg_free,lower.tail=FALSE)
  if (Qbvg >= CHIbvg){
    print("H2:FAILURE")
  } else {
    print("H2:PASS")
  }
  print(deg_free)
  # Third test: same center but different covariance structure
  test_stats = test.cov(n1,n2,m1,m2,s1,s2)
  # Will print the number of iterations needed to calculate covariance matrices in the console
  h1_cov = matrix(test_stats[[1]],ncol=length(m1)) # Converts the list output into usable square matrices
  h2_cov = matrix(test_stats[[2]],ncol=length(m1))
  Qcvg = n1*log(det(h1_cov)/det(s1))+n2*log(det(h2_cov)/det(s2))
  deg_free = D-1
  CHIcvg = qchisq(0.67,deg_free,lower.tail=FALSE)
  if (Qcvg >= CHIcvg){
    print("H3:FAILURE")
  } else {
    print("H3:PASS")
  }
  print(deg_free)
  results = array(0,dim=c(3,2)) # Initialize an array for storage
  results[,1] = c(Qavg,Qbvg,Qcvg)
  results[,2] = c(CHIavg,CHIbvg,CHIcvg)
  return(results)
}