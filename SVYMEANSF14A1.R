

rm(list=ls())
library(stochvol)
library(coda)
set.seed(0) # for replication
options(digits=17)
vt   = read.table('vytF14A1.txt',sep = '\t')
T    = dim(vt)[1]
N    = dim(vt)[2]
for (i in 1:N){
  if(min(log(vt[,i]^2))== -Inf){
    vt[,i] = vt[,i] + 0.00001 #offset to avoid taking log of zero
  }
}

# Run MCMC algorithm and store draws
S    = 50000
burn = 50000
m    = matrix(0,T+3,N)
g    = matrix(0,4,N)
for (i in 1:N){
  draws  = svsample(vt[,i],draws=S,burnin=burn,quiet=TRUE,thinpara=10,thinlatent=10)
  all1  = data.frame(draws$para[[1]])
  all2 = data.frame(draws$latent[[1]])
  all1 = colMeans(all1)
  all2 = t(all2)
  all2 = rowMeans(all2)
  all2 = unname(all2)
  all1 = unname(all1)
  
  m[1:3,i] = all1[1:3]
  for (j in 1:618){
    m[j+3,i] = all2[j]
  }
  
  g[,i]  = geweke.diag(draws$para[[1]][,-4])$z
  name   = sprintf('svydraws%d.txt',i)
  #	write(t(all),file=name,ncolumn=dim(all)[2])
}
out = rbind(m,g) #include Geweke statistics
write(t(out),file='svymeansF14A1.txt',ncolumn=dim(out)[2])
