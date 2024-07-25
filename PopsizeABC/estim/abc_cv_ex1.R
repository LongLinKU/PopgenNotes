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

source(paste0(importDir,'/generations.R'))
source(paste0(importDir,'/myfunctions.R'))


nval=100 # number of tested histories

args = commandArgs(trailingOnly=TRUE)
print(args)
infile_params = fs::path_abs(args[1])
infile_stat = fs::path_abs(args[2])
pop=args[3]
n=as.numeric(args[4])
mac = ceiling(n*0.1)
tol = as.numeric(args[5]) # 0.1
outputDir = args[6]


# loads simulated samples, all SNP for AFS stats and MAF 20% for LD stats
source(paste0(importDir,'/load_simu.R'))

# choice of summary stats
nb_afs=n/2-mac+2
ind_afs=1:nb_afs
ind_ld=nb_afs+(1:(nb_dist-1)) # the LD statistic corresponding to the shortest distance is removed
ind_ibs=nb_afs+nb_dist+(1:(nb_m*nb_prob))
ind_stat=c(ind_afs,ind_ld) # only AFS and LD statistics

# abc
res=cv4abc(param=log10_pop(params)[,-1],sumstat=stat[,ind_stat],nval=nval,tols=c(tol),statistic="median",method="neuralnet") #a tolerance of order 0.001 would give much better results, but this would require more simulated samples.
par.estim=(res$estim)[[1]]
par.true=res$true
errors=colSums((par.estim-par.true)**2)/(length(par.true[,1])*var(par.true[,dim(par.true)[2]]))

# plot
pdf(paste0(outputDir,'/ex4_',pop,'_n',n,'_mac_',mac,'_tols_',tol,'.cv.pdf'),height=4,width=4)
par(mar=c(3,3,1,2),cex=0.7,mgp=c(1.5,0.5,0))
plot(NA,xlim=c(0,5.5),ylim=c(0,1),xlab="generations before present (log scale)",ylab="population size prediction error ",xaxt="n")
axis(1,at=0:5,labels=c("0","10","100","1,000","10,000","100,000"))
y=errors[-1]
points(c(generations,5.5),c(y,y[nb_times]),type="s",lwd=2)
dev.off()


