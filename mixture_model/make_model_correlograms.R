# make model correlograms for each file in the directory
# created 2019-01-29

aalist = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

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





#