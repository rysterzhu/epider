#******************************************************************************
#	Filename:	go_analysis_barplot.r
#	Author:DingGui Tao
#******************************************************************************
setwd('~//')
tt <- read.delim('.txt', header=T)
tt <- tt[1:20,]
pval <- format(tt$PValue, digits=2, scientific=T)
term <- unlist(strsplit(as.character(tt$Term),'~'))[seq(2,40,2)]
par(oma=c(15,2,2,2))
b <- barplot(tt$Count) # col=heat.colors(20)
text(b[,1],par('usr')[3]-0.025,adj=1,srt=45, labels=term, xpd=NA)
text(b[,1],tt$Count/2,pval,srt=90)


#******************************************************************************
#	Filename:	GO_analysis_barplot.r
#	Author:	Leah Sun
#	Date:	2016-05-19
#******************************************************************************
library(Himsc)
library(devEMF)
setwd('~//')
tt <- read.delim('.txt', header=T)
term <- unlist(strsplit(as.character(tt$Term), '~'))[seq(2,22,2)]
emf('GO_anlysis_result.emf',width=8, height=7)
bp <- barplot(-log10(PValue), tt$Count,horiz=TRUE, col=c('lightblue'), xlim=c(0,4),
              xlab='-Log10(P-Value)', space=0.25, cex.axis=1.3, cex.lab=2, lwd=2)
text(0, bp, term, pos=4, cex=1.3)
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()

################################################################################
setwd("C:/Users/cate/Desktop/graduate/")
tiff(file='loss_GO.tiff', height=3200, width=3200, units = "px", res = 600, compression='lzw')

#enhancer_go=matrix(c('nucleosome assembly',42,'8.83E-21','dendrite morphogenesis',21,'2.9e-4','nervous system development',18,'1.3e-5',
#                    'sensory perception of pain',31,'5.8e-3','axon guidance',30,'1.4e-10','proteolysis',27,'8.5e-3')
#                   ,ncol=3,byrow = T)

tt <-  read.csv('loss_go.csv', header=FALSE)
enhancer_go <- as.matrix(tt)
enhancer_go[,2:3]=as.numeric(enhancer_go[,2:3])
bp=barplot(-log10(as.numeric(enhancer_go[,3])),horiz = T,col='#56B4E9',xlab='-log10(PValue)')
text(x = 0,y=bp,labels = rev(enhancer_go[,1]),pos = 4)
axis(2,at=bp,labels = as.numeric(enhancer_go[,2]),las=2,tick=F,line=-0.9)
#legend('bottomright',legend = c('neuron','glia'),fill=c('lightblue','pink'))

dev.off()
