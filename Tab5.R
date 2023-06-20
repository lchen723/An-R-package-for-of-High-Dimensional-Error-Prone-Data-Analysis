
install.packages("SIMEXBoost")
library(SIMEXBoost)
R=c(0.1,0.3,0.5)
data = read.table("bankruptcy_data.csv",sep=",",head=T)
data = data[,-94]
Y = data[,95]
Y = as.numeric(Y)
Xstar = t(as.matrix(data[,-95]))
Xstar=scale(Xstar)

p = dim(Xstar)[1]
n = dim(Xstar)[2]

set.seed(202270)
##naive
EST<-NULL
naive<-Boost_VSE(Y,t(Xstar),type="binary",Iter = 50,Lambda = 0)$BetaHat
EST<-cbind(EST,naive)
##linear
for(i in R){
  correctL <- SIMEXBoost(Y,t(Xstar),zeta=c(0,0.25,0.5,0.75,1),B=5, type="binary",
                         sigmae=diag(i,p), Iter=50, Lambda=0,Extrapolation="linear")$BetaHatCorrect
  
  EST = rbind(EST,t(correctL))
}
##Quadratic

for(i in R){
  correctQ <- SIMEXBoost(Y,t(Xstar),zeta=c(0,0.25,0.5,0.75,1),B=5, type="binary",
                         sigmae=diag(i,p), Iter=50, Lambda=0,Extrapolation="quadratic")$BetaHatCorrect

  EST = rbind(EST,t(correctQ))
}

data<-as.data.frame(t(EST))
data<-na.omit(data)
naive<-data[,1]
beta0.1L<-data[,2]
beta0.3L<-data[,3]
beta0.5L<-data[,4]
beta0.1Q<-data[,5]
beta0.3Q<-data[,6]
beta0.5Q<-data[,7]

output<-function(X){
  b<-data.frame(cbind(which(X !=0),X[which(X!=0)]))
  return(b)
}
resultnaive<-output(naive)
result0.1L<-output(beta0.1L)
result0.3L<-output(beta0.3L)
result0.5L<-output(beta0.5L)
result0.1Q<-output(beta0.1Q)
result0.3Q<-output(beta0.3Q)
result0.5Q<-output(beta0.5Q)





