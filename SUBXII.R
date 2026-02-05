#-------------------------------------------------
#Unit Sine Burr XII
#SUBXII(mu,sigma)
#-------------------------------------------------
#Log Likelihood Function
LSUBURR<-expression(log(-(pi*mu*sigma*(((-log(y))^sigma+1)^mu+2)*((-log(y))^sigma+1)^(-2*mu-1)*
                            (-log(y))^sigma*cos((pi*(1/((-log(y))^sigma+1)^mu+1))/
                                                  (4*((-log(y))^sigma+1)^mu)))/(4*y*log(y))))
#-------------------------------------------------
#Score Vector
dfdm<-D(LSUBURR,"mu")
d2fd2m<-D(dfdm,"mu")
dfds<-D(LSUBURR,"sigma")
d2fd2s<-D(dfds,"sigma")
d2fdmds<-D(dfdm,"sigma")
#-------------------------------------------------
#GAMLSS Function
SUBXII<-function(mu.link="log",sigma.link="log") 
{
  mstats <- checklink("mu.link", "SUBXII", substitute(mu.link), 
                      c("log","identity"))
  dstats <- checklink("sigma.link", "SUBXII", substitute(sigma.link), 
                      c("log","identity"))
  structure(list(family = c("SUBXII","Unit Sine Burr XII"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
                 nopar = 2, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 dldm = function(y, mu, sigma) {
                   dldm <- eval(dfdm)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   d2ldm2 <- eval(d2fd2m)
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {
                   dldd <- eval(dfds)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   d2ldd2 = eval(d2fd2s)
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   d2ldmdd = eval(d2fdmds)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 #--------------------------------------------------------------
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dSUBXII(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pSUBXII",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(0.5, length(y))),
                 mu.valid = function(mu) all(mu > 0), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ), 
  class = c("gamlss.family", "family"))
}
#-------------------------------------------------
#Probability density function
dSUBXII<-function(y,mu,sigma,log=FALSE)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  fy1<--(pi*mu*sigma*(((-log(y))^sigma+1)^mu+2)*((-log(y))^sigma+1)^(-2*mu-1)*(-log(y))^sigma*cos((pi*(1/((-log(y))^sigma+1)^mu+1))/(4*((-log(y))^sigma+1)^mu)))/(4*y*log(y))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  return(fy)
}
#-------------------------------------------------
#Cumulative distribution function
pSUBXII<-function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("q must be between 0 and 1", "\n", ""))
  G<-(1+(-log(q))^sigma)^(-mu)
  cdf1<-sin((pi/4)*G*(G+1))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  return(cdf)
}
#-------------------------------------------------
#Quantile function
qSUBXII<-function(p,mu,sigma)
{
  term <- ((-1 + sqrt(1 + (16 * asin(p)) / pi)) / 2)^(-1 / mu) - 1
  y<-exp(-term^(1 / sigma))
  return(y)
}
#-------------------------------------------------
#SUBXII Random Number Generation Function
rSUBXII<-function(n,mu,sigma)
{
  u<-runif(n)
  term <- ((-1 + sqrt(1 + (16 * asin(u)) / pi)) / 2)^(-1 / mu) - 1
  y<-exp(-term^(1 / sigma))
  return(y)
}
#-------------------------------------------------
#GSUBXII Random Number Generation Function
r_gsubxii<-function(n,mu=NULL,sigma=NULL,beta1=NULL,beta2=NULL,X1=NULL,X2=NULL)
{
  glink<-make.link("log")
  u<-runif(n)
  if(!is.null(beta1) & !is.null(X1)){
    mu<-glink$linkinv(X1%*%beta1)
  }
  if(!is.null(beta2) & !is.null(X2)){
    sigma<-glink$linkinv(X2%*%beta2)
  }
  term <- ((-1 + sqrt(1 + (16 * asin(u)) / pi)) / 2)^(-1 / mu) - 1
  y<-exp(-term^(1 / sigma))
  return(y)
}
#-------------------------------------------------
#Examples using the gamlss package
library(gamlss)
#-------------------------------------------------
#SUBXII Distribution
y<-rSUBXII(n=100,mu=6,sigma=2)
df<-as.data.frame(y)
plot.ts(df$y)
hist(df$y)
mod<-gamlss::gamlss(y~1,sigma.formula=~1,family=SUBXII(mu.link="identity",sigma.link = "identity"),data=df)
summary(mod)
#-------------------------------------------------
#GSUBXII Model
#True Parameters
n<-100
beta1_true<-c(1.5,-2.2,1.8)
beta2_true<-c(0.8,1.8)

#Covariáveis
k1<-length(beta1_true)-1
X1<-cbind(rep(1,n),matrix(runif(n*k1),nrow=n,ncol=k1))
k2<-length(beta2_true)-1
X2<-cbind(rep(1,n),matrix(runif(n*k2),nrow=n,ncol=k2))
y<-r_gsubxii(n=n,mu=mu_true,sigma=sigma_true,beta1=beta1_true,beta2=beta2_true,X1=X1,X2=X2)
plot.ts(y)
hist(y)

df<-data.frame(y=y,X1=X1[,2:3],X2=X2[,2])
mod<-gamlss::gamlss(y~X1.1+X1.2,sigma.formula=~X2,family=SUBXII(),data=df)
summary(mod)
#-------------------------------------------------

