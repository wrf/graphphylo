# created 2018-10-22 by WRF

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/supplements/feuda2017/whelan2017/filtered_Day_CAT_GTR2_RFd.txt"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

chainrfd = read.table(inputfile,header=FALSE,sep="")

pdf(file=outputfile, width=6, height=5)
par( mar=c(4.5,4.5,2,1) )
plot(chainrfd[,1], chainrfd[,3], type='l', col="#01665e", lwd=2, xlab="Iteration", ylab="Robinson-Foulds distance", frame.plot=FALSE, main=inputfile, cex.lab=1.4, cex.axis=1.4)

dev.off()


#