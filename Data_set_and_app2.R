#Sivini,2006
dados<-c(8,4.67,5,1.17,5.50,1.67,1.25,1.17,1.17,1.33,2.17,1.50,2.33,1.33,1.17,6,1.17,6)


hist(dados)
require(AdequacyModel)
TTT(dados)
library(MASS)
require(GenSA)
require(survival)
require(moments)

n=length(dados)
n
is.numeric(dados)

truehist(dados,
         ylim=c(0,0.6),
         col = "white",ylab="f(x)",xlab = "x",nbins = 6)

summary(dados)
tabela=round(cbind(mean(dados), median(dados), sd(dados), var(dados),skewness(dados), kurtosis(dados), min(dados), max(dados)),4)
tabela

fit.sa2<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2]))) #minus the loglik  
  lower <- c(0,0) 
  upper <- c(10,10)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=7000,max.time=3))
  return(out[c("value","par","counts")])
}


## 
fit.sa3<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3]))) #minus the loglik  
  lower <- c(0,0,0) 
  upper <- c(100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=10000,max.time=3))
  return(out[c("value","par","counts")])
}

fit.sa4<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3],x[4]))) #minus the loglik  
  lower <- c(0,0,0,0) 
  upper <- c(100,100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}

##################
#WFGa
##################
pdf.WNFGa = function(x,b,alpha,beta){
  H<-pgamma(x,alpha,beta)
  h<-dgamma(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1730)
fit.sa3(dados,pdf.WNFGa)

pdf_WNFGa = function(par, x){
  b = par[1]
  alpha = par[2]
  beta = par[3]
  
  H<-pgamma(x,alpha,beta)
  h<-dgamma(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_WNFGa = function(par, x){
  b = par[1]
  alpha = par[2]
  beta = par[3]
  
  H<-pgamma(x,alpha,beta)
  h<-dgamma(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fa=1-exp(-a*(G/(1-G))^(b))
  fa
}

set.seed(1730)
resultsWNFGa  = goodness.fit(pdf = pdf_WNFGa , cdf = cdf_WNFGa ,
                             starts = c( 0.0535549, 37.1446638, 12.2865824),
                             data = dados, method = "BFGS", domain = c(0, Inf),
                             mle = NULL);resultsWNFGa





####################################
#EGa
###############################

pdf.EGa = function(x,a,alpha,beta){
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  fd=a*g*G^(a-1)
  fd
}
set.seed(1730)
fit.sa3(dados,pdf.EGa)

pdf_EGa = function(par, x){
  a = par[1]
  alpha = par[2]
  beta = par[3]
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  fd=a*g*G^(a-1)
  fd
}
cdf_EGa = function(par, x){
  a = par[1]
  alpha = par[2]
  beta = par[3]
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  fa=G^a
  fa
}

set.seed(1730)
resultsEGa  = goodness.fit(pdf = pdf_EGa , cdf = cdf_EGa ,
                           starts = c(100,   0.03031326,   0.42346729),
                           data = dados, method = "BFGS", domain = c(0, Inf),
                           mle = NULL);resultsEGa
#
chutes<-c(100,   0.03031326,   0.42346729)
fit.gtnhEGa<- fitdistr(dados,pdf.EGa,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-8,1e-8,1e-8),maxit=10000))
fit.gtnhEGa


###########################
#MO-Ga
###########################
pdf.moGa<-function(x,a,alpha,beta){
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  res<-a*g/(1-(1-a)*(1-G))^2
  return(res)
}
set.seed(1730)
fit.sa3(dados,pdf.moGa)

pdf_moGa<-function(par,x){
  a<-par[1]
  alpha<-par[2]
  beta<-par[3]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  
  res<-a*g/(1-(1-a)*(1-G))^2
  return(res)
}

#acumulada
cdf_moGa<-function(par,x){
  a<-par[1]
  alpha<-par[2]
  beta<-par[3]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  
  
  res<-1 - a*(1-G)/(1-(1-a)*(1-G))
  return(res)
}
set.seed(1730)
resultsmoGa  = goodness.fit(pdf = pdf_moGa , cdf = cdf_moGa,
                            starts = c(0.1407559, 2.9452956, 0.5663591),
                            data = dados, method = "BFGS", domain = c(0, Inf),
                            mle = NULL);resultsmoGa

chutes<-c(0.1407559, 2.9452956, 0.5663591)
fit.gtnhmoGa<- fitdistr(dados,pdf.moGa,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-8,1e-10,1e-6),maxit=10000))
fit.gtnhmoGa



#####################
#GGa
############
pdf.gga = function(x,a,alpha,beta){
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}

set.seed(1729)
fit.sa3(dados,pdf.gga)

pdf_gga = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}


cdf_gga = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  
  fa=pgamma(-log(1-G),a)
  
}

set.seed(1729)
resultsgga  = goodness.fit(pdf = pdf_gga , cdf = cdf_gga ,
                           starts = c(1.4577696, 1.4867262, 0.7043512),
                           data = dados, method = "BFGS", domain = c(0, Inf),
                           mle = NULL);resultsgga

chutes<-c(1.4577696, 1.4867262, 0.7043512)
fit.gtnhgga<- fitdistr(dados,pdf.gga,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-6,1e-8,1e-6),maxit=10000))
fit.gtnhgga


########################################################
#Kw-Ga
#######################################################

pdfkwga<-function(x,a,b,alpha,beta){
  
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
}
set.seed(1729)
fit.sa4(dados,pdfkwga)

pdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
}
cdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwga = goodness.fit(pdf = pdf_Kwga, cdf = cdf_Kwga,
                           starts = c( 12.3409824,  0.1232018,  2.1147009,  4.6014399), 
                           data = dados, method = "BFGS", domain = c(0, Inf),
                           mle = NULL);resultsKwga

chutes<-c(12.3409824,  0.1232018,  2.1147009,  4.6014399)
fit.gtnhkwga<- fitdistr(dados,pdfkwga,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-10,1e-10,1e-6,1e-12),maxit=10000))
fit.gtnhkwga



#################
# Beta Ga
#################
bga.pdf=function(x,a,b,alpha,beta){
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
set.seed(1729)
fit.sa4(dados,bga.pdf)

cdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  pbeta(G, a, b)
}
pdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  G<-pgamma(x,alpha,beta)
  g<-dgamma(x,alpha,beta)
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbga = goodness.fit(pdf = pdfbga, cdf = cdfbga,
                          starts = c( 0.05364416,  0.07140238, 53.82061592, 17.13915499), 
                          data = dados, method = "BFGS", domain = c(0, Inf),
                          mle = NULL);resultsbga





########
#Ga
#######
pdfga<-function(x,alpha,beta){
  
  
  g<-dgamma(x,alpha,beta)
  
}
set.seed(1729)
fit.sa2(dados,pdfga)

pdf_ga = function(par, x){
  
  alpha = par[1]
  beta = par[2]
  
  g<-dgamma(x,alpha,beta)
  
}
cdf_ga = function(par, x){
  
  alpha = par[1]
  beta = par[2]
  G<-pgamma(x,alpha,beta)
  
  
}
set.seed(1729)
resultsga = goodness.fit(pdf = pdf_ga, cdf = cdf_ga,
                         starts = c( 2.0799395, 0.7117664), 
                         data = dados, method = "BFGS", domain = c(0, Inf),
                         mle = NULL);resultsga


#GRAFICOS
#densities
truehist(dados,
         ylim=c(0,0.6),
         col = "white",ylab="f(x)",xlab = "x",nbins = 6)


curve(pdf_WNFGa(resultsWNFGa$mle,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(pdf_moGa(fit.gtnhmoGa$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(pdf_Kwga(fit.gtnhkwga$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=3)
curve(pdfbga(resultsbga$mle,x),add=TRUE, lwd = 3, col="green",lty=4)

legend(4,.5, legend = c("WFGa", "MOGa","KwGa", "BGa"),
       col = c("red","blue","yellow","green"), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")

####cdfs
km<- survfit(Surv(dados) ~ 1) #Kaplan-Meier
plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l")#cdf empririca
curve(cdf_WNFGa(resultsWNFGa$mle,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(cdf_moGa(fit.gtnhmoGa$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(cdf_Kwga(fit.gtnhkwga$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=3)
curve(cdfbga(resultsbga$mle,x),add=TRUE, lwd = 3, col="green",lty=4)

legend(5,0.6, legend = c("WFGa", "MOGa","KwGa", "BGa"),
       col = c("red","blue","yellow","green"), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")

