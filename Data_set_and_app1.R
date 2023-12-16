#Óbitos mulheres idade fértil segundo Região - http://tabnet.datasus.gov.br/cgi/tabcgi.exe?sim/cnv/mat10uf.def
#Período: 2016-2020
#valores divididos por 100

dados<-c(68.44,210.02,313.33,91.09,57.77,57.41,179.63,
         270.82,85.20,49.52,57.62,179.68,264.30,85.08,50.25,
         56.79,183.94,266.58,87.62,48.73,55.02,189.17,282.18,
         92.16,52.94)

n=length(dados)
n

hist(dados)
require(AdequacyModel)
TTT(dados, col="red", lwd=2.5, grid=TRUE, lty=2)
library(MASS)
require(GenSA)
require(survival)
require(moments)

n=length(dados)
n
is.numeric(dados)


summary(dados)
tabela=round(cbind(mean(dados), median(dados), sd(dados), var(dados),skewness(dados), kurtosis(dados), min(dados), max(dados)),4)
tabela

fit.sa2<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2]))) #minus the loglik  
  lower <- c(0,0) 
  upper <- c(1000,100000)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,maxit=7000,max.time=3))
  return(out[c("value","par","counts")])
}


##
fit.sa3<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3]))) #minus the loglik  
  lower <- c(0,0,0) 
  upper <- c(100,1000,100)
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

################################
#BP
#########################

require(extraDistr)
bp.pdf=function(x,a,b){
  g<-dbetapr(x,a,b)
  
}

set.seed(1729)
fit.sa2(dados,bp.pdf)

pdf_bp = function(par, x){
  a = par[1]
  b = par[2]
  g<-dbetapr(x,a,b)
  
}

cdf_bp = function(par, x){
  a = par[1]
  b = par[2]
  
  G<-pbetapr(x,a,b)
  
}
set.seed(1729)
resultsbp = goodness.fit(pdf = pdf_bp, cdf = cdf_bp,
                         starts = c(4.452255 ,  2.759698), 
                         data = dados, method = "C", domain = c(0, Inf),
                         mle = NULL);resultsbp




############################

#WFBP

############################
require(extraDistr)
pdf.WNFBP = function(x,b,alpha,beta){
  h=dbetapr(x,alpha,beta)
  H=pbetapr(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
set.seed(1730)
fit.sa3(dados,pdf.WNFBP)



pdf_WNFBP = function(par, x){
  b = par[1]
  alpha = par[2]
  beta = par[3]
  
  h=dbetapr(x,alpha,beta)
  H=pbetapr(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fd=a*b*g*(G^(b-1)/(1-G)^(b+1))*exp(-a*(G/(1-G))^(b))
  fd
}
cdf_WNFBP = function(par, x){
  b = par[1]
  alpha = par[2]
  beta = par[3]
  
  h=dbetapr(x,alpha,beta)
  H=pbetapr(x,alpha,beta)
  
  g <- h*(1-H)^H*(H/(1-H)-log(1-H))
  G <- 1-(1-H)^H
  
  a<-1
  fa=1-exp(-a*(G/(1-G))^(b))
  fa
}

set.seed(1730)
resultsWNFBP  = goodness.fit(pdf = pdf_WNFBP , cdf = cdf_WNFBP ,
                             starts = c(0.1905959, 1000,  9.7030240),
                             data = dados, method = "C", domain = c(0, Inf),
                             mle = NULL);resultsWNFBP


chutes<-c(0.1905959, 1000,  9.7030240)
fit.gtnhWNFBP<- fitdistr(dados,pdf.WNFBP,start=list(b=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-6,1e-10,1e-10),maxit=10000))
fit.gtnhWNFBP





####################################
#EBP
###############################

pdf.EBP = function(x,a,alpha,beta){
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  
  fd=a*g*G^(a-1)
  fd
}
set.seed(1730)
fit.sa3(dados,pdf.EBP)

pdf_EBP = function(par, x){
  a = par[1]
  alpha = par[2]
  beta = par[3]
  
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  
  fd=a*g*G^(a-1)
  fd
}
cdf_EBP = function(par, x){
  a = par[1]
  alpha = par[2]
  beta = par[3]
  
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  
  fa=G^a
  fa
}

set.seed(1730)
resultsEBP  = goodness.fit(pdf = pdf_EBP , cdf = cdf_EBP ,
                           starts = c(0.1904711, 1000, 4.6555988),
                           data = dados, method = "C", domain = c(0, Inf),
                           mle = NULL);resultsEBP

chutes<-c(0.1904711, 1000, 4.6555988)
fit.gtnhEBP<- fitdistr(dados,pdf.EBP,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-10,1e-10,1e-12),maxit=10000))
fit.gtnhEBP




###########################
#MO-BP
###########################
pdf.mobp<-function(x,a,alpha,beta){
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  res<-a*g/(1-(1-a)*(1-G))^2
  return(res)
}
set.seed(1730)
fit.sa3(dados,pdf.mobp)

pdf_mobp<-function(par,x){
  a<-par[1]
  alpha<-par[2]
  beta<-par[3]
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  
  res<-a*g/(1-(1-a)*(1-G))^2
  return(res)
}


cdf_mobp<-function(par,x){
  a<-par[1]
  alpha<-par[2]
  beta<-par[3]
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  
  
  res<-1 - a*(1-G)/(1-(1-a)*(1-G))
  return(res)
}
set.seed(1730)
resultsmobp  = goodness.fit(pdf = pdf_mobp , cdf = cdf_mobp,
                            starts = c(0.5709063, 2.6135960,   2.4761439),
                            data = dados, method = "C", domain = c(0, Inf),
                            mle = NULL);resultsmobp

chutes<-c(0.5709063, 2.6135960,   2.4761439)
fit.gtnhmobp<- fitdistr(dados,pdf.mobp,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-8,1e-8,1e-8),maxit=10000))
fit.gtnhmobp



#####################
#GBP
############
pdf.gBP = function(x,a,alpha,beta){
  
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}

set.seed(1729)
fit.sa3(dados,pdf.gBP)

pdf_gBP = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  
  g=dbetapr(x,alpha,beta)
  G=pbetapr(x,alpha,beta)
  
  fd=g/gamma(a)*(-log(1-G))^(a-1)
  fd
}


cdf_gBP = function(par, x){
  
  a= par[1]
  alpha = par[2]
  beta = par[3]
  
  g=(beta*alpha^beta)/(x^(beta+1))
  G<-1-(alpha/x)^(beta)
  fa=pgamma(-log(1-G),a)
  
}

set.seed(1729)
resultsgBP  = goodness.fit(pdf = pdf_gBP , cdf = cdf_gBP ,
                           starts = c(2 ,12  , 3.320332),
                           data = dados, method = "C", domain = c(0, Inf),
                           mle = NULL);resultsgBP

chutes<-c(2 ,12  , 3.320332)
fit.gtnhgBP<- fitdistr(dados,pdf.gBP,start=list(a=chutes[1],alpha=chutes[2],beta=chutes[3]),control=list(ndeps=c(1e-12,1e-12,1e-6),maxit=10000))
fit.gtnhgBP




########################################################
#Kw-bp
#######################################################
kwbp.pdf=function(x,alpha,beta,a,b){
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  fd<-alpha*beta*g*G^(alpha-1)*(1-G^alpha)^(beta-1)
  fd
}
set.seed(1729)
fit.sa4(dados,kwbp.pdf)

pdf_Kwbp = function(par, x){
  alpha = par[1]
  beta = par[2]
  a = par[3]
  b = par[4]
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  alpha*beta*g*G^(alpha-1)*(1-G^alpha)^(beta-1)
}

cdf_Kwbp = function(par, x){
  alpha = par[1]
  beta = par[2]
  a = par[3]
  b = par[4]
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  fa<-1-(1-G^alpha)^beta
}
set.seed(1729)
resultsKwbp = goodness.fit(pdf = pdf_Kwbp, cdf = cdf_Kwbp,
                           starts = c(1,0.2485105, 1, 6.3183283), 
                           data = dados, method = "C", domain = c(0, Inf),
                           mle = NULL);resultsKwbp


#################
# Beta bp
#################
bbp.pdf=function(x,alpha,beta,a,b){
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  fd=(1/beta(alpha,beta))*g*(G^(alpha-1))*(1-G)^(beta-1) 
  fd
}

set.seed(1729)
fit.sa4(dados,bbp.pdf)

cdfbbp=function(par,x)
{
  alpha = par[1]
  beta = par[2]
  a = par[3]
  b = par[4]
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  pbeta(G, alpha, beta)
}
pdfbbp=function(par,x)
{
  alpha = par[1]
  beta = par[2]
  a = par[3]
  b = par[4]
  g<-dbetapr(x,a,b)
  G<-pbetapr(x,a,b)
  (1/beta(alpha,beta))*g*(G^(alpha-1))*(1-G)^(beta-1)
}
set.seed(1729)
resultsbbp = goodness.fit(pdf = pdfbbp, cdf = cdfbbp,
                          starts = c(10,  0.2335827, 100,  6.6500398),
                          data = dados, method = "C", domain = c(0, Inf),
                          mle = NULL);resultsbbp

############
chutes<-c(10,  0.2335827, 100,  6.6500398)
fit.gtnhbbp<- fitdistr(dados,bbp.pdf,start=list(alpha=chutes[1],beta=chutes[2],a=chutes[3], b=chutes[4]),control=list(ndeps=c(1e-8,1e-10,1e-8,1e-8),maxit=10000))
fit.gtnhbbp



####################################################

#graficos
truehist(dados,
         ylim=c(0,.022),
         col = "white",ylab="f(x)",xlab = "x",nbins = 10)
curve(pdf_WNFBP(fit.gtnhWNFBP$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(pdf_EBP(fit.gtnhEBP$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(pdf_gBP(fit.gtnhgBP$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=3)
curve(pdfbbp(fit.gtnhbbp$estimate,x),add=TRUE, lwd = 3, col="green",lty=4)

legend(170,.015, legend = c("WFBP", "EBP", "GBP", "BBP"),
       col = c("red","blue","yellow","green"), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")

km<- survfit(Surv(dados) ~ 1) #Kaplan-Meier
plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l")#cdf empririca
curve(cdf_WNFBP(fit.gtnhWNFBP$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(cdf_EBP(fit.gtnhEBP$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(cdf_gBP(fit.gtnhgBP$estimate,x),add=TRUE, lwd = 3, col="yellow",lty=3)
curve(cdfbbp(fit.gtnhbbp$estimate,x),add=TRUE, lwd = 3, col="green",lty=4)
legend(200,0.6, legend = c("Empirical", "WFBP", "EBP", "GBP", "BBP"),
       col = c("black","red","blue","yellow","green"), lwd = 3 ,lty = c(1,1,2,3,4) , bty ="n")


