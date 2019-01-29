# plot direct comparison of buried and exposed substitution matrices for EX2 model
# model by Le, Lartillot, and Gascuel 2008
# script created 2019-01-29

colorset = c( colorRampPalette(c("#543005ee","#ffffd488"),alpha=TRUE)(9),colorRampPalette(c("#ddd4f088","#253494ee"),alpha=TRUE)(8) )
#plot(1:length(colorset),1:length(colorset), pch=16, cex=5, col=colorset)
basecolors = c("#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788")

bur_dat = read.table("~/git/graphphylo/mixture_model/bur_EX2_model.hpi.txt", sep=" ", fill=TRUE, col.names=aalist[s_hyd$ix])
#bur_dat
#sum(rowSums(bur_dat,na.rm=TRUE))
exp_dat = read.table("~/git/graphphylo/mixture_model/exp_EX2_model.hpi.txt", sep=" ", fill=TRUE, col.names=aalist[s_hyd$ix])
#exp_dat

pdf(file="~/git/graphphylo/mixture_model/EX2_model_bur-exp_difference.pdf", width=7, height=7)
par(mar=c(4,0.5,1,3))
plot(0,0,type="n",xlim=c(0,20),ylim=c(1,20),axes=FALSE, xlab="", ylab="")
axis(1,at=1:20,labels=FALSE, cex.axis=1.4)
mtext(aalist[s_hyd$ix], side=1, at=1:20, cex=1.4, line=1)
axis(4,at=1:20,labels=FALSE, cex.axis=1.4)
mtext(aalist[s_hyd$ix], side=4, at=1:20, cex=1.4, line=1)
# plot difference of each row
for (i in 1:19) {
	i
	as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])
	abs(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i]))
	log(abs(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])))
	4+log(abs(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])))/2
	pointscaling = 4+log( abs(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])))
	#floor(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])+9)
	colorindex = floor(as.numeric(bur_dat[i,1:i]-exp_dat[i,1:i])+9)
	colorindex[colorindex<1]=1
	colorindex[colorindex>17]=17
	pointcol = colorset[ colorindex ]
	points(rep(i+1,i), 1:i, cex=abs(pointscaling)/1.5, pch=16, col=pointcol)
}
pointscaling = as.numeric(bur_dat[20,1:20] - exp_dat[20,1:20])*20
as.numeric(bur_dat[20,1:20] - exp_dat[20,1:20])*20
colorindex = floor( as.numeric(bur_dat[20,1:20] - exp_dat[20,1:20]) * 20 + 4)
points(1:20, 1:20, cex=abs(pointscaling), pch=15, lwd=2, col=basecolors[colorindex])
# make legend
points( rep(1.5,7), seq(14,17,0.5), cex=c(3.5,2,1,1,1,2,3.5), pch=16, col=colorset[c(2,4,6,8,12,14,16)])
#legend(1,17, legend=c("Transition preferred\nwhen buried","","","","","","Transition preferred\nwhen exposed"), pt.cex=c(3,2,1,0.5,1,2,3), pch=16, col=colorset[c(2,4,6,9,12,14,16)], bty='n')
#rect( rep(1,17), seq(13,17,0.25), rep(2,17), seq(13.25,17.25,0.25), col=colorset)
text(2,17, "Transition preferred\nwhen buried", cex=1.5, pos=4)
text(2,14, "Transition preferred\nwhen exposed", cex=1.5, pos=4)
legend(1,11.5,legend=c("Common exposed", "Common buried"), pch=15, pt.cex=c(2,2), col=basecolors[c(1,6)], cex=1.3, bty='n')
# add model name
text(1,20, "EX2 Model: Buried - Exposed", cex=2, pos=4)
dev.off()


