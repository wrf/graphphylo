# make model correlograms for each file in the directory
# created 2019-01-29

aalist = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V") # alpha by AA name

#aalist = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") # alpha by letter
#aalist = c("D", "P", "E", "N", "K", "R", "Q", "S", "G", "H", "T", "A", "C", "Y", "M", "V", "W", "L", "I", "F") # hpi

# light teal to blue scheme for overall transition frequency
colorset = c("#7fcdbb99", "#41b6c4aa", "#1d91c0aa", "#225ea8cc", "#081d58ee")
#plot(1:length(colorset),1:length(colorset), pch=16, cex=5, col=colorset)

# purple color range for displaying overall base frequency 
basecolors = c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#756bb1", "#54278f")

# data from https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
# which itself derives from:
# Relationship of sidechain hydrophobicity and a‐helical propensity on the stability of the single‐stranded amphipathic a‐helix
# https://doi.org/10.1002/psc.310010507
hydrophobicity = c(41,-14,-28,-55,49,-10,-31,0,8,99,97,-23,74,100,-46,-5,13,96,63,76)
# sort and return index, which will reorder from hydrophilic (Asp) to most hydrophobic (Phe)
s_hyd = sort(hydrophobicity, index.return=TRUE)

modeldir = "~/git/graphphylo/mixture_model/"
filelist = dir(modeldir,pattern="*.hpi.txt")
for (modelfile in filelist) {
	joinfilename = paste0(modeldir, modelfile)
	#print(joinfilename)
	modeldat = read.table(joinfilename, sep=" ", fill=TRUE, col.names=aalist[s_hyd$ix])
	outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",joinfilename,perl=TRUE)
	print(outputfile)
	pdf(file=outputfile, width=7, height=7)
	par(mar=c(4,0.5,1,3))
	plot(0,0,type="n",xlim=c(0,20),ylim=c(1,20),axes=FALSE, xlab="", ylab="")
	axis(1,at=1:20,labels=FALSE, cex.axis=1.4)
	mtext(aalist[s_hyd$ix], side=1, at=1:20, cex=1.4, line=1)
	axis(4,at=1:20,labels=FALSE, cex.axis=1.4)
	mtext(aalist[s_hyd$ix], side=4, at=1:20, cex=1.4, line=1)
	for (i in 1:19) {
		pointscaling = (6 + log(as.numeric(modeldat[i,1:i])))/2
		colorindex = floor(pointscaling)+1
		colorindex[colorindex<1]=1
		points(rep(i+1,i), 1:i, cex=pointscaling, pch=16, col=colorset[colorindex])
	}
	pointscaling = (6 + log(as.numeric(modeldat[20,1:20])))/1.3
	colorindex = floor(pointscaling)+1
	points(1:20,1:20, cex=pointscaling, col=basecolors[colorindex], pch=15)
	text(1,20, basename(modelfile), cex=2, pos=4)
	legend(2,14,legend=c("AA rare", "AA common"), pch=15, pt.cex=c(1,3), col=basecolors[c(2,4)], cex=1.3)
	legend(2,18,legend=c("Transition infrequent", "Transition frequent"), pch=16, pt.cex=c(1,3), col=colorset[c(1,4)], cex=1.3)
	dev.off()
}



# default order LG model
colorset = colorRampPalette(c("#7fcdbb99", "#225ea8cc", "#020928ee"))(6)
joinfilename = "~/git/graphphylo/mixture_model/LG_model.txt"
modeldat = read.table(joinfilename, sep=" ", fill=TRUE, col.names=aalist, row.names=c(aalist[2:20],aalist[1]) )
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",joinfilename,perl=TRUE)
print(outputfile)
pdf(file=outputfile, width=7, height=7)
par(mar=c(4,3,1,0.5))
plot(0,0,type="n",xlim=c(0,20),ylim=c(1,20),axes=FALSE, xlab="", ylab="")
axis(1,at=1:20,labels=FALSE, cex.axis=1.4)
mtext(aalist, side=1, at=1:20, cex=1.4, line=1)
axis(2,at=1:20,labels=FALSE, cex.axis=1.4, line=-1)
mtext(rev(aalist), side=2, at=1:20, cex=1.4, line=0)
for (i in 1:19) {
	pointscaling = (6 + log(as.numeric(modeldat[i,1:i])))/2
	colorindex = floor(pointscaling)+1
	colorindex[colorindex<1]=1
	points(1:i, rep(20-i,i), cex=pointscaling, pch=16, col=colorset[colorindex])
}
pointscaling = (6 + log(as.numeric(modeldat[20,1:20])))/1.3
colorindex = floor(pointscaling)+1
points(1:20,20:1, cex=pointscaling, col=basecolors[colorindex], pch=15)
text(20,20, basename(joinfilename), cex=2, pos=2)
legend(20,14,legend=c("AA rare", "AA common"), pch=15, pt.cex=c(1,3), col=basecolors[c(2,4)], cex=1.3, xjust=1)
legend(20,18,legend=c("Transition infrequent", "Transition frequent"), pch=16, pt.cex=c(1,3), col=colorset[c(1,4)], cex=1.3, xjust=1)
dev.off()


# default order WAG model
# point scaling parameters vary substantially, range appears to be from appx 1 to 999
colorset = colorRampPalette(c("#7fcdbb99", "#225ea8cc", "#020928ee"))(8)
joinfilename = "~/git/graphphylo/mixture_model/WAG_model.txt"
modeldat = read.table(joinfilename, sep=" ", fill=TRUE, col.names=aalist, row.names=c(aalist[2:20],aalist[1]) )
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",joinfilename,perl=TRUE)
print(outputfile)
pdf(file=outputfile, width=7, height=7)
par(mar=c(4,3,1,0.5))
plot(0,0,type="n",xlim=c(0,20),ylim=c(1,20),axes=FALSE, xlab="", ylab="")
axis(1,at=1:20,labels=FALSE, cex.axis=1.4)
mtext(aalist, side=1, at=1:20, cex=1.4, line=1)
axis(2,at=1:20,labels=FALSE, cex.axis=1.4, line=-1)
mtext(rev(aalist), side=2, at=1:20, cex=1.4, line=0)
for (i in 1:19) {
	pointscaling = log(as.numeric(modeldat[i,1:i]),base=4)
	colorindex = floor(pointscaling)+1
	colorindex[colorindex<1]=1
	points(1:i, rep(20-i,i), cex=pointscaling, pch=16, col=colorset[colorindex])
}
pointscaling = (6 + log(as.numeric(modeldat[20,1:20])))/1.3
colorindex = floor(pointscaling)+1
points(1:20,20:1, cex=pointscaling, col=basecolors[colorindex], pch=15)
text(20,20, basename(joinfilename), cex=2, pos=2)
legend(20,14,legend=c("AA rare", "AA common"), pch=15, pt.cex=c(1,3), col=basecolors[c(2,4)], cex=1.3, xjust=1)
legend(20,18,legend=c("Transition infrequent", "Transition frequent"), pch=16, pt.cex=c(1,3), col=colorset[c(1,4)], cex=1.3, xjust=1)
dev.off()











#