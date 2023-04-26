# Calculate the common mean from two given data populations
common.mean <- function(n1,n2,s1,s2,m1,m2){
  # n = number of observations per sample [scalar]
  # s = sample covariance matrix [D-1,D-1]
  # m = sample mean [1,D-1]
  a = solve((n1*solve(s1))+(n2*solve(s2)))
  b = (n1*solve(s1)%*%m1)+(n2*solve(s2)%*%m2)
  c_mean = a%*%b
  return(t(c_mean))
}