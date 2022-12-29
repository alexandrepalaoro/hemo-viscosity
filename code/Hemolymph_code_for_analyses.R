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


## Removing rustic for the sensitivity analysis that comes later on! ###

no.rustic = hawkm[!(row.names(hawkm) %in% ("Manduca rustica")),]

no.rus.tree = keep.tip(pruned.h,rownames(no.rustic))

plot(no.rus.tree)
dev.off()

### All good to proceed with the analyses ###

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

#### NOW, WE RUN THE SENSITIVITY ANALYSES ####
## The idea is to remove Manduca rustica and see how much it is influencing the pattern we found with the full dataset ##

no.rusOU = phylolm(m.visc ~ forewing.length + I(forewing.length^2) , data=no.rustic , phy = no.rus.tree, "OUfixedRoot",
                      upper.bound = 3, measurement_error = T)

no.rusLB = phylolm(m.visc ~ forewing.length + I(forewing.length^2) , data=no.rustic , phy = no.rus.tree, "lambda",
                      lower.bound = 1e-10)

no.rusBM = phylolm(m.visc ~ forewing.length + I(forewing.length^2), data=no.rustic, phy = no.rus.tree, "BM", 
                      measurement_error = T)

no.rusOU$aic; no.rusLB$aic; no.rusBM$aic

## Again, BM is the best fit model.

no.rus.step = phylostep(m.visc ~ forewing.length + I(forewing.length^2), data=no.rustic, phy = no.rus.tree, model = "BM", 
                      measurement_error = T, boot = 10^4, direction = 'both')

summary(no.rus.step)

#### SENSITVITY ANALYSIS FIGURE ####
png("figures/viscosity-sensitvt.png",res=600,units='mm',h=180,w=220)

plot(m.visc~forewing.length,data=no.rustic,las=1,bty='l',
     pch=21,cex=2.7,ylim=c(1,2.8),
     bg=c("grey50","orange")[as.numeric(as.factor(no.rustic$group))],
     ylab = "Mean viscosity (mPa*s)", xlab= "Mean forewing length (mm)")

arrows(x0=no.rustic$forewing.length, x1=no.rustic$forewing.length,
       y0=no.rustic$m.visc-no.rustic$visc.std, y1=no.rustic$m.visc+no.rustic$visc.std, 
       code=3, angle=90, length=0.1, lwd = 2)

points( x = no.rustic$forewing.length, y = no.rustic$m.visc, cex = 2.8,
        pch = 21, bg = c("grey50","orange")[as.numeric(as.factor(no.rustic$group))])

curve(no.rus.step$coefficients[1]+(no.rus.step$coefficients[2]*x),add=T,lwd=5,col=alpha('grey30',0.5),lty=2)

text(m.visc~forewing.length, data = no.rustic, labels = label, font = 2, pos = 4, cex = 0.8, offset = 0.6)

dev.off()

##### DATA ON THE DISTANCE BETWEEN MUSCLE FIBERS ####

#pdf("figures/viscosity-slits.pdf",w=12,h=8)
plot(mean.slit ~ forewing.length, data = hawkm, las = 1, bty = 'l',
     cex = 2, pch = 21, bg = c("grey","orange")[as.numeric(as.factor(hawkm$group))],
     ylab = "Average distance between longitudinal muscle fibers (mm)",
     xlab = "Forewing length (cm)", ylim = c(0,0.5))

arrows(x0=hawkm$forewing.length, x1=hawkm$forewing.length,
       y0=hawkm$mean.slit-hawkm$sd.slit, y1=hawkm$mean.slit+hawkm$sd.slit, 
       code=3, angle=90, length=0.1, lwd = 2)

points( x = hawkm$forewing.length, y = hawkm$mean.slit, cex = 2.3,
        pch = 21, bg = c("grey50","orange")[as.numeric(as.factor(hawkm$group))])


text(mean.slit ~ forewing.length, data = hawkm, labels = label, font = 2, pos = 4.5, cex = 0.8)
abline(v = 4, lwd = 5, lty = 2, col = 'grey')
#dev.off()

## NOTE: I added the corresponding labels in photoshop. 
