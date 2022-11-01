library(ape)
library(phylolm)
library(phytools)
library(geiger)
library(treeio)
library(RColorBrewer)
library(scales)


##### LOADING DATASET #####

hawkm = read.csv("data/hemolymph_data.csv",h=T,sep=';')

## TREE LOADING, PRUNING, AND ADDING ####

phy = read.beast("phylo/Hawkmoth.202tx.BEAST.tre")

pruned.h = keep.tip(phy@phylo,hawkm$sp.in.tree)

hawkm$species = gsub("_", " ", hawkm$species)

pruned.h$tip.label

pruned.h$tip.label = c( "Agrius cingulata",  "Ceratomia catalpae", "Dolba hyloeus" , 
                        "Eumorpha pandorus", "Enyo lugubris", "Hemaris diffinis", 
                        "Hyles lineata", "Manduca quinquemaculata", "Manduca rustica" , 
                        "Manduca sexta", "Paratrea plebeja", "Xylophanes tersa")


rownames(hawkm) = hawkm$species

name.check(pruned.h,hawkm)

### This the tree we used for analyses
#png("figures/phylo1.png",res=600,w=250,h=250,units='mm')
plot.phylo(ladderize(pruned.h),label.offset = 2, edge.width = 2);axisPhylo()
mtext("Time (my)", side = 1, line = 3, at = 21.35)
#dev.off()



#### PGLSs #####

### BODY VOLUME ANALYSES ####

### First, we built the full models with different evolutionary models, and compared their AIC

vol.highOU = phylolm(m.visc ~ volume + I(volume^2) , data=hawkm , phy = pruned.h, model = "OUfixedRoot",
                     upper.bound = 4, measurement_error = T)

vol.highLB = phylolm(m.visc ~ volume + I(volume^2) , data=hawkm , phy = pruned.h, model = "lambda",
                     lower.bound = 1e-10)

vol.highBM = phylolm(m.visc ~ volume + I(volume^2), data=hawkm, phy = pruned.h, model = "BM", 
                     measurement_error = T)


vol.highOU$aic; vol.highLB$aic; vol.highBM$aic

### Lambda and BM have the lowest AIC value. Thus, we keep the BM model because it has fewer parameters 
### when compared to the lambda.

## Now, stepwise model selection with bootstrap to calculate the CIs.

vol.step = phylostep(m.visc ~ volume + I(volume^2), data=hawkm, phy = pruned.h, model = "BM", 
                     measurement_error = T, boot = 10^4, direction = 'both')

summary(vol.step)


### FOREWING LENGTH ANALYSES ####

wing.highOU = phylolm(m.visc ~ forewing.length + I(forewing.length^2) , data=hawkm , phy = pruned.h, "OUfixedRoot",
                      upper.bound = 3, measurement_error = T)

wing.highLB = phylolm(m.visc ~ forewing.length + I(forewing.length^2) , data=hawkm , phy = pruned.h, "lambda",
                      lower.bound = 1e-10)

wing.highBM = phylolm(m.visc ~ forewing.length + I(forewing.length^2), data=hawkm, phy = pruned.h, "BM", 
                      measurement_error = T)

wing.highOU$aic; wing.highLB$aic; wing.highBM$aic

## Again, BM is the best fit model.

wing.step = phylostep(m.visc ~ forewing.length + I(forewing.length^2), data=hawkm, phy = pruned.h, model = "BM", 
                      measurement_error = T, boot = 10^4, direction = 'both')

summary(wing.step)


### PREPARING THE PHENOGRAM AND THE ANCESTRAL RECONSTRUCTION ####

m.visc = setNames(hawkm$m.visc, rownames(hawkm))

plot(pruned.h)
nodelabels()

pruned2 = paintSubTree(pruned.h, node=14, state='2')

cols<-c("grey80","orange"); names(cols)<-1:2

phenogram(pruned2,m.visc,colors=cols,ylab="Mean viscosity (mPa*s)",xlab="Time since root (my)",ftype='i', 
          ylim = c(1.4,2.3))


visc.anc = contMap(ladderize(pruned.h),m.visc,res=1000,plot=F)

visc.anc = setMap(visc.anc, colors=brewer.pal(9, name = "Greens"))

#### FIGURE #####

#png("figures/viscosity-figs4.png",res=600,units='mm',h=380,w=380)
par(mfrow=c(2,2),cex=1.1)

plot(m.visc~volume,data=hawkm,las=1,bty='l',
     pch=21,cex=3,bg=c("grey50","orange")[as.numeric(as.factor(hawkm$group))],
     ylim = c(1, 2.8), ylab = "Mean viscosity (mPa*s)", xlab= "Mean body volume (cm³)")

arrows(x0=hawkm$volume, x1=hawkm$volume,
       y0=hawkm$m.visc-hawkm$visc.std, y1=hawkm$m.visc+hawkm$visc.std, 
       code=3, angle=90, length=0.1, lwd = 2)

points( x = hawkm$volume, y = hawkm$m.visc, cex = 3,
        pch = 21, bg = c("grey50","orange")[as.numeric(as.factor(hawkm$group))])

legend("topleft","a)",bty='n',cex=1.2)


plot(m.visc~forewing.length,data=hawkm,las=1,bty='l',
     pch=21,cex=3,ylim=c(1,2.8),
     bg=c("grey50","orange")[as.numeric(as.factor(hawkm$group))],
     ylab = "Mean viscosity (mPa*s)", xlab= "Mean forewing length (mm)")

arrows(x0=hawkm$forewing.length, x1=hawkm$forewing.length,
       y0=hawkm$m.visc-hawkm$visc.std, y1=hawkm$m.visc+hawkm$visc.std, 
       code=3, angle=90, length=0.1, lwd = 2)

points( x = hawkm$forewing.length, y = hawkm$m.visc, cex = 3,
        pch = 21, bg = c("grey50","orange")[as.numeric(as.factor(hawkm$group))])

curve(wing.step$coefficients[1]+(wing.step$coefficients[2]*x)+ (wing.step$coefficients[3]*x^2),add=T,lwd=5,col=alpha('grey30',0.5),lty=2)
legend("topleft","b)",bty='n',cex=1.2)

phenogram(pruned2,colors=cols,m.visc,ylab="Mean viscosity (mPa*s)",xlab="Time since root (my)",
          ftype='i',ylim=c(1.3,2.3), las = 1)
legend("topleft","c)",bty='n',cex=1.2)

plot(visc.anc, legend = 28, leg.txt = "Mean viscosity (mPa*s)")

legend("topleft","d)",bty='n',cex=1.2)

#dev.off()
