library(evd)


alpha<-2
beta<-0.5
gamma=beta/alpha
#C=1/gamma(1-beta)
t=seq(1:10000)
U=runif(t)
W=rstable(t,alpha=beta,beta=0,gamma=1,delta=0)
#W=(U/C)^(-1/beta)-(beta/(beta-1))*C^(1/beta)
Z=rfrechet(t,1,0,alpha)
J<-W^{1/alpha}*Z
S=cumsum(W)
plot(S,J,type="h")

estimates <- function(TT, idxJ, KK){
  # Creates a dataframe with log-moment estimates
  # of the tail and scale parameters
  # TT:   arrival times of the magnitudes JJ
  # idxJ: index vector which orders JJ decreasingly
  # KK:   a vector of k's for which the estimates are to be calculated.
  #       The threshold is set at the k'th largest observation.
  returned_df <- ldply(.data = KK, function(k){
    WW <- get_durations(TT, idxJ, k)
    est <- ml.par.est.delta(WW,0.05)                                          #in der AppVersion ist hier noch ein zus??tzlicher Parameter??
    row <- c(k,est$nu, est$CInu, est$delta, est$CIdelta)
    return(row)
  })
  names(returned_df) <- c("k", "tail", "tailLo", "tailHi",
                          "scale", "scaleLo", "scaleHi")
  return(returned_df)
}

#uncoupedFall stability Plots
Ju=rfrechet(t,0,1,alpha)
df<-data.frame(times=W,magnitudes=Ju)
df$idxJ <- order(df$magnitudes, decreasing = TRUE)
estims<-estimates(W,df$idxJ,seq(150,10000,1))
delta0 <- estims$scale * ((estims$k))^(beta)
par(mfrow=c(2,1))
plot(estims$k,estims$tail,type="l")
plot(estims$k, delta0, type='l', main = "ML scale")


#Laplace transform inversion

fun<-function(s)
{
  exp(-s^0.5)
}

x <- seq(0.01,100, 0.25)
sinvals <- vapply(x, iv.opC, complex(1), L.FUN = fun)
plot(x, Re(sinvals), type = "l")
lines(x, dstable(x, alpha = 0.5, beta = 0),col="blue")


plot(Re(sinvals),dstable(x, alpha = 0.5, beta = 0))