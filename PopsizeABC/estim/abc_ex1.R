rm(list=ls())
library(abc)
whereDir <- function(){
    # function to get directory where scripts are, so accompanying functions can be sourced even if the script is run from outside
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}
importDir = whereDir()

args = commandArgs(trailingOnly=TRUE)
print(args)

gen_time = as.numeric(args[8])
# loads simulated samples, MAF 20% for both AFS and LD stats

source(paste0(importDir,'/generations.R'))
source(paste0(importDir,'/myfunctions.R'))


infile_params = fs::path_abs(args[1])
infile_stat = fs::path_abs(args[2])
popMy=args[3]
n=as.numeric(args[4])
mac = ceiling(n*0.1)
source(paste0(importDir,"/load_simu.R"))

# loads observed samples
infile_obs=fs::path_abs(args[5])
source(paste0(importDir,"/load_obs.R"))

tol = as.numeric(args[6]) # 0.1
outputDir = args[7]


# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
ind_ibs=nb_afs+nb_dist+(1:(nb_m*nb_prob))
ind_stat=c(ind_afs,ind_ld) # only AFS and LD statistics

# abc
abc_res=abc(obs[ind_stat],log10_pop(params)[,-1],stat[,ind_stat],tol=tol,method="neuralnet") # a tolerance of order 0.001 would give much better results, but this would require more simulated samples.
abc_estim=summary(abc_res,print=FALSE)[3,] # median
#abc_estim=summary(abc_res,print=FALSE)[5,] # mode

# plot
pdf(paste0(outputDir,'/ex3_',popMy,'_n',n,'_mac_',mac,'_tols_',tol,'.pdf'),height=4,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5),ylim=c(2,5),xlab="years before present (log scale)",ylab="effective population size (log scale)",axes=F)
axis(1,at=0:5,labels=c("0","10","100","1,000","10,000","100,000"))
axis(2,at=c(2,log10(300),3,log10(3000),4,log10(30000),5),labels=c("100","300","1,000","3,000","10,000","30,000","100,000"))
lines(years,abc_estim[-1],type="s",lwd=2) # the first term of the vector is removed because it corresponds to the recombination rate
abc_q=summary(abc_res,print=FALSE,intvl=0.9)[2,] # 5% quantile
print(abc_q)
lines(years,abc_q[-1],type="s",lwd=2,lty=3)
abc_q=summary(abc_res,print=FALSE,intvl=0.9)[6,] # 95% quantile
print(abc_q)
lines(years,abc_q[-1],type="s",lwd=2,lty=3)
dev.off()

# plot detailed information about the estimation of all parameters
detailedFile=paste0(outputDir,'/ex3_',popMy,'_n',n,'_mac_',mac,'_tols_',tol,'.est.details.pdf')

plot(abc_res,param=log10_pop(params)[,-1],file='abc_details',subsample=20)
file.rename('abc_details.pdf',detailedFile)
write.table(data.frame(years=years,abc_estim=abc_estim[-1]),paste0(outputDir,'/ex3_',popMy,'_n',n,'_mac_',mac,'_tols_',tol,'.est.median.txt'),sep="\t", quote=FALSE, col.names = T, row.names = F) 




