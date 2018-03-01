# plot phylobayes trace log likelihoods

# usage:

# Rscript plot_phylobayes_traces.R hif-pb1.trace hif-pb2.trace

# add optional max iterations next:

# Rscript plot_phylobayes_traces.R hif-pb1.trace hif-pb2.trace 10000

args = commandArgs(trailingOnly=TRUE)
trace1file = args[1]
trace2file = args[2]
othermax = args[3]
optionname = args[4]

#trace1file = "~/data/project/supplements/feuda2017/whelanD20/CAT_GTR1.trace"
#trace2file = "~/data/project/supplements/feuda2017/whelanD20/CAT_GTR2.trace"
trace1dat = read.table(trace1file, header=FALSE, sep="\t")
trace2dat = read.table(trace2file, header=FALSE, sep="\t")

#trace1file = "~/data/project/supplements/feuda2017/whelan2017/Day_CAT_GTR1.trace"
#trace2file = "~/data/project/supplements/feuda2017/whelan2017/Day_CAT_GTR2.trace"

#trace1file = "~/est/feuda2017/bzxdp-feuda_et_al_2017-e970b278280d/CHANG/CHANG_ANALYSES_SUMMARY_FILES/LG1.trace"
#trace2file = "~/est/feuda2017/bzxdp-feuda_et_al_2017-e970b278280d/CHANG/CHANG_ANALYSES_SUMMARY_FILES/LG2.trace"

### CHECK FOR ALTERNATE HEADER LINE
if ( trace1dat[1,1]=="iter" ) {
trace1dat = read.table(trace1file, header=TRUE, sep="\t")
trace2dat = read.table(trace2file, header=TRUE, sep="\t")
}


paste("trace 1 of length", length(trace1dat[,1]) )
paste("trace 2 of length", length(trace2dat[,1]) )

logln1 = trace1dat[,4]
minmaxlogln1 = range(logln1)
minmaxlogln1
firstsat1 = minmaxlogln1[1] - (minmaxlogln1[1] - minmaxlogln1[2]) * 0.99
paste("first value at 99% of max is", firstsat1 )

tracemax = max( c(length(trace1dat[,1]),length(trace2dat[,1])) )
maxlength = min( min( which(logln1 >= firstsat1))+min(5000,round(length(trace1dat[,1])/10)), tracemax )
#maxlength = 500
### MAKE SURE MAX CANNOT GO MORE THAN TRACE LENGTH
if (!is.na(othermax)) {
	paste("given number of iterations is", othermax )
	maxlength = min( c(as.numeric(othermax), tracemax ) )
}
paste("using", maxlength, "as max iterations in plot")

llrange  = range( pretty( range(trace1dat[,4])) )
serange  = range( pretty( range( c( trace1dat[,8],trace2dat[,8]) ) ) )
lenrange = range( pretty( range( c( trace1dat[,5],trace2dat[,5]) ) ) )
alpharange = range( pretty( range( c( trace1dat[,6],trace2dat[,6]) ) ) )
noderange = range( pretty( range( c( trace1dat[,7],trace2dat[,7]) ) ) )
statalpha = range( pretty( range( c( trace1dat[,9],trace2dat[,9]) ) ) )

### LG AND WAG MODELS DO NOT HAVE 10 COLUMNS
if (dim(trace1dat)[2] < 10) {
rrrange  = statalpha # meaning use alpha instead
} else {
rrrange  = range( pretty( range( c( trace1dat[,10],trace2dat[,10]) ) ) )
}

legendlabels = c( paste("Trace 1 (",length(trace1dat[,1]),"iter.)"), paste("Trace 2 (",length(trace2dat[,1]),"iter.)") )

### MAKE GRAPH
outputfile = paste(trace1file,".pdf",sep="")
#pdf(file="~/data/project/supplements/feuda2017/whelan2017/Day_CAT_GTR1.trace.pdf", width=8, height=7)
pdf(file=outputfile, width=12, height=9)
par(mfrow=c(2,3), mar=c(4.5,4.5,4,2) )
plot( trace1dat[,1], trace1dat[,4], type='l', xlim=c(0, maxlength), ylim=llrange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="log(L)" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,4]+10, lwd=3, col="#0571b099"  )
#text(maxlength*0.55, llrange[1] + diff(llrange) * 0.6, paste("Final iterations =",tracemax) , cex=1.5)
legend(maxlength*0.21, llrange[1] + diff(llrange) * 0.4, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

plot( trace1dat[,1], trace1dat[,8], type='l', xlim=c(0, maxlength), ylim=serange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="statent" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,8], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, serange[1] + diff(serange) * 0.6, paste("Final iterations =",tracemax) , cex=1.5)
legend(maxlength*0.21, serange[1] + diff(serange) * 0.9, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

### CHECK FOR OPTIONAL LABEL NAME
if (!is.na(optionname)) {
#mtext( paste(optionname, trace1file,"+",trace2file), side=3, cex=1.4, at=0-maxlength*0.2, line=1)
mtext( paste(optionname, trace1file,"+",trace2file), side=3, cex=1.4, line=1)
} else {
#mtext( paste(trace1file,"+",trace2file), side=3, cex=1.4, at=0-maxlength*0.2, line=1)
mtext( paste(trace1file,"+",trace2file), side=3, cex=1.4, line=1)
}

plot( trace1dat[,1], trace1dat[,5], type='l', xlim=c(0, maxlength), ylim=lenrange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="length" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,5], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, lenrange[1] + diff(lenrange) * 0.6, paste("Final iterations =",max( c(length(trace1dat[,1]),length(trace2dat[,1])) )) , cex=1.5)
legend(maxlength*0.21, lenrange[1] + diff(lenrange) * 0.4, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

plot( trace1dat[,1], trace1dat[,6], type='l', xlim=c(0, maxlength), ylim=alpharange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="alpha" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,6], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, alpharange[1] + diff(alpharange) * 0.6, paste("Final iterations =",max( c(length(trace1dat[,1]),length(trace2dat[,1])) )) , cex=1.5)
legend(maxlength*0.21, alpharange[1] + diff(alpharange) * 0.4, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

plot( trace1dat[,1], trace1dat[,7], type='l', xlim=c(0, maxlength), ylim=noderange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="Nmode" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,7], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, noderange[1] + diff(noderange) * 0.6, paste("Final iterations =",max( c(length(trace1dat[,1]),length(trace2dat[,1])) )) , cex=1.5)
legend(maxlength*0.21, noderange[1] + diff(noderange) * 0.4, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

mtext( paste("Final iterations =",tracemax), side=3, cex=1.4, line=1)

### MAKE DIFFERENT PLOTS FOR RRENT AND ALPHA
if (dim(trace1dat)[2] < 10) {
plot( trace1dat[,1], trace1dat[,9], type='l', xlim=c(0, maxlength), ylim=statalpha, lwd=3, col="#f1a34099", xlab="Iterations", ylab="statalpha" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,9], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, statalpha[1] + diff(statalpha) * 0.6, paste("Final iterations =",max( c(length(trace1dat[,1]),length(trace2dat[,1])) )) , cex=1.5)
legend(maxlength*0.21, statalpha[1] + diff(statalpha) * 0.4, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)
} else {
plot( trace1dat[,1], trace1dat[,10], type='l', xlim=c(0, maxlength), ylim=rrrange, lwd=3, col="#f1a34099", xlab="Iterations", ylab="rrent" , main="", cex.lab=1.4 )
lines( trace2dat[,1], trace2dat[,10], lwd=3, col="#0571b099"  )
#text(maxlength*0.75, rrrange[1] + diff(rrrange) * 0.6, paste("Final iterations =",max( c(length(trace1dat[,1]),length(trace2dat[,1])) )) , cex=1.5)
legend(maxlength*0.21, rrrange[1] + diff(rrrange) * 0.9, legend=legendlabels, col=c("#f1a340", "#0571b0"), lwd=2, cex=1.4)

}


dev.off()




#