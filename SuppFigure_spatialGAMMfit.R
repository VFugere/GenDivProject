## Vincent Fugere 2019
## Supp Figure showing GAMM fit for all taxa

colz <- c(1,'#E69F00','#56B4E9','#009E73')

pdf('~/Desktop/FigS1.pdf',width=8,height=10,pointsize = 8,onefile = T)

par(mfrow=c(4,4),cex=1)

plot_smooth(m1, view="D", cond=list('lat' = 0, 'long' = 0), col=colz[1], rm.ranef=T, hide.label = T,xlab='mean spatial distance (D)',ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
legend('topleft',bty='n',legend=c('scale = 10 km','100 km','1 000 km','10 000 km'),pch=16,col=colz)

plot_smooth(m2, view="D", cond=list('lat' = 0, 'long' = 0), col=colz[2], rm.ranef=T, hide.label = T,xlab='mean spatial distance (D)',ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m3, view="D", cond=list('lat' = 0, 'long' = 0), col=colz[3], rm.ranef=T, hide.label = T,xlab='mean spatial distance (D)',ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m4, view="D", cond=list('lat' = 0, 'long' = 0), col=colz[4], rm.ranef=T, hide.label = T,xlab='mean spatial distance (D)',ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))

plot_smooth(m1, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=colz[1], rm.ranef=T, hide.label = T,xlab=expression(latitude~(absolute)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m2, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=colz[2], rm.ranef=T, hide.label = T,xlab=expression(latitude~(absolute)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m3, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=colz[3], rm.ranef=T, hide.label = T,xlab=expression(latitude~(absolute)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m4, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=colz[4], rm.ranef=T, hide.label = T,xlab=expression(latitude~(absolute)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))

plot_smooth(m1, view="hd", cond=list('lat' = 0, 'long' = 0), col=colz[1], rm.ranef=T, hide.label = T,xlab=expression(human~density~(people~km^-2)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m2, view="hd", cond=list('lat' = 0, 'long' = 0), col=colz[2], rm.ranef=T, hide.label = T,xlab=expression(human~density~(people~km^-2)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m3, view="hd", cond=list('lat' = 0, 'long' = 0), col=colz[3], rm.ranef=T, hide.label = T,xlab=expression(human~density~(people~km^-2)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m4, view="hd", cond=list('lat' = 0, 'long' = 0), col=colz[4], rm.ranef=T, hide.label = T,xlab=expression(human~density~(people~km^-2)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))

plot_smooth(m1, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=colz[1], rm.ranef=T, hide.label = T,xlab=expression(land~use~intensity),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m2, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=colz[2], rm.ranef=T, hide.label = T,xlab=expression(land~use~intensity),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m3, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=colz[3], rm.ranef=T, hide.label = T,xlab=expression(land~use~intensity),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))
plot_smooth(m4, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=colz[4], rm.ranef=T, hide.label = T,xlab=expression(land~use~intensity),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0,0.008))

dev.off()