
#######################################
#######################################
# ANCESTRAL RANGE RECONSTRUCTION - RASE
#######################################
#######################################
library(rase)
library(raster)
library(ape)

setwd("C:/Users/carol/Documents/Workshop_meninas/biogeobears-meninas/files_originais")

#transform shape in owin format
results<-readRDS("tree_mapa.rds")
tree.sigmodontinae<-results[[1]]
sigmodontinae<-results[[2]]
sigmodontinae2<-raster::aggregate(sigmodontinae,by="binomial")


owin.shapes.regular<-shape.to.rase(sigmodontinae2) # regular projection
owin.shapes<-owin.shapes.regular


#y<- as(sigmodontinae2, "SpatialPolygons")
#p <- slot(y, "polygons")
#v <- lapply(p, function(z) { SpatialPolygons(list(z)) })
#str(v[[1]])
#winlist <- lapply(v, as.owin)

# assign names to polys & include underline
sig.names<-unique(sigmodontinae2$binomial)
sig.names1<-sub(" ","_",sig.names,fixed=TRUE) 
sig.names1 == tree.sigmodontinae$tip.label
sig.names1<-sort(sig.names1)

tree.sigmodontinae$tip.label
length(sig.names1)
#
winlist<-name.poly(owin.shapes,tree.sigmodontinae,poly.names = sig.names1)

# run rase
rase_results<-rase(tree.sigmodontinae,winlist,niter=12000)




library(polyCub)
str(rase_results)
fix(rase_results)
#write.table(rase_results,"rase_results.txt")
#rase_results<-read.table("rase_results.txt")


# Use the amazing 'coda' package to explore the MCMC
require(coda)
# post-MCMC handling
rasemcmc <- coda::mcmc(rase_results)
str(rasemcmc)
summary_values_rase<-summary(rasemcmc)
str(summary_values_rase)

summary_values_rase$quantiles
nrow(summary_values_rase$statistics)

#options(max.print=1000)
summary_values_rase$statistics[517:520,]
#plot the traces for all the parameters 
plot(rasemcmc)
plot(rasemcmc[,519:520],trace=FALSE,density=TRUE,smooth=TRUE)
densplot(rasemcmc)




#### RENAN a partir daqui é pra ser os plots, pelo que entendi!




#plot
# transform in 3D
#df3<-data.for.3d(rase_results,tree.sigmodontinae,winlist) # rase or mcmc
# plot 3D
#require(rpanel)
#phylo.3d(df3,z.scale=1)
# add the polygons representing the tip distributions
#add.polygons(df3)
# add the posterior density at each node of the 3d tree
#add.dens(df3, rase_results, nlevels = 10,z.scale=1)


#Slices # corta no valor desejado, ex: 12.65Ma e nós "vivos" ali
slice_results<-rase.slice(tree.sigmodontinae,slice=12.65,res=rase_results,winlist,niter=100)
slice_results<-rase.slice(tree.sigmodontinae,slice=10,res=rasemcmc,winlist,niter=100)
slice_results<-rase.slice(tree.sigmodontinae,slice=8,res=rasemcmc,winlist,niter=100)
### os dois próximos dão erro!!! by Carol
slice_results<-rase.slice(tree.sigmodontinae,slice=6,res=rasemcmc,winlist,niter=100) ### erro
slice_results<-rase.slice(tree.sigmodontinae,slice=4,res=rasemcmc,winlist,niter=100) ### erro

slice_results
rsl=slice_results
str(slice_results)
slicemcmc <- coda::mcmc(slice_results)
slicemcmc

#df_slice<-data.for.3d(slice_results,tree.sigmodontinae,winlist) # rase or mcmc
#phylo.3d(df_slice,z.scale=5)
#add.polygons(df_slice)
#add.dens(df_slice, slice_results, z.scale = 5)


#############################################################
#############################################################
#############################################################
###
#Mapping Slices - From Quintero
###
library(rgeos)
library(maptools)
library(raster)
library(ks)
library(RColorBrewer)

# tree to see dates and nodes

plot(tree.sigmodontinae,cex=0.3,type="fan")
nodelabels(cex=0.3, frame = "circle") # nó 374 é Akodontini
axisPhylo()

# check who is alive (second column is less than `number of species+1`)
tree.slice(tree.sigmodontinae, 12.65) # [,1] node ancestral; [,2]node descendant
tree.slice(tree.sigmodontinae, 10) 

plot(tree.sigmodontinae)
axisPhylo()
nodelabels(cex=0.3, frame = "circle")
edgelabels(cex=0.3, frame = "circle")

#color palettes in order of polygons
col.pals = c('Oranges', 'BuPu', 'PuRd', 'Greys', 'YlGn','YlGnBu','YlOrBr', 'Greens', 'Reds', 'Blues', 'BuGn','Oranges', 'BuPu', 'PuRd', 'Greys', 'YlGn','YlGnBu','YlOrBr', 'Greens', 'Reds', 'Blues', 'BuGn','Oranges', 'BuPu', 'PuRd', 'Greys', 'YlGn','YlGnBu','YlOrBr', 'Greens', 'Reds', 'Blues', 'BuGn','Oranges', 'BuPu', 'PuRd', 'Greys', 'YlGn','YlGnBu','YlOrBr', 'Greens', 'Reds', 'Blues', 'BuGn','Oranges', 'BuPu', 'PuRd', 'Greys', 'YlGn','YlGnBu','YlOrBr', 'Greens', 'Reds', 'Blues', 'BuGn')

#plot(shape.neot.beh)



plot(shape.neot)
library(makeFlow)
library(scales)

# Loop for plotting each polygon
for (i in 1:(ncol(rsl)/2)) {
  
  # make `x` & `y` object
  df = data.frame(rsl[,i], rsl[,i+(ncol(rsl)/2)])
  
  # estimate Highest Posterior Interval
  hh = Hpi(df, binned = TRUE)*1
  
  # make 2D dimensional smoothing
  dd = kde(df, H = hh)
  
  # create contour lines for polygons
  cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                    z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                  length = 20))
  
  # set color

  #col1 = addAlpha(brewer.pal(3, col.pals[i]), 0.4)
  #col1 = brewer.pal(3, col.pals[i])
  col1 <- alpha(brewer.pal(3, col.pals[i]), 0.4)
    #col1=c("pink","blue","yellow")

  # plot three stacked polygons for each level of the HPI
  polygon(cc[[16]]$x, cc[[16]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
  polygon(cc[[17]]$x, cc[[17]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
  polygon(cc[[18]]$x, cc[[18]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
}

####
library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(letsR)
shape.neot <- readOGR(choose.files()) # Folder: Dados Espaciais > Book Chapter Neotropical Diversification
shape.neot

#plot(shape.neot.beh)
plot(shape.neot)

# Loop for plotting the one polygon
for (i in 1:(ncol(rsl)/2)) {
  
  # make `x` & `y` object
  df = data.frame(rsl[,2], rsl[,2+(ncol(rsl)/2)])
  
  # estimate Highest Posterior Interval
  hh = Hpi(df, binned = TRUE)*1
  
  # make 2D dimensional smoothing
  dd = kde(df, H = hh)
  
  # create contour lines for polygons
  cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                    z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                  length = 20))
  
  # set color
  col1 = alpha(brewer.pal(3, col.pals[i]), 0.4)
  #col1=c("pink","blue","yellow")
  
  # plot three stacked polygons for each level of the HPI
  polygon(cc[[14]]$x, cc[[14]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
  polygon(cc[[16]]$x, cc[[16]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
  polygon(cc[[18]]$x, cc[[18]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
}




#############################################################
#############################################################
# Direto com rase_results (1 reconstrução para cada nó)
rase_results # duas últimas colunas são valores de sigma (-1)
rase_results[,2] 
# Ex:
plot(tree.sigmodontinae)
nodelabels(cex=0.3) # nó 261 é o ancestral (raiz), nó 262 é Sigmodontini
# plot map
#plot(shape.neot.beh)
plot(shape.neot)
# Loop for plotting the one polygon
for (i in 1:(ncol(rsl)/2)) {
  
  # subset of rase results
  rase_results_subset<-rase_results[1:10,]
  
  # make `x` & `y` object
  df = data.frame(rase_results_subset[,261], rase_results_subset[,261+(ncol(rase_results)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
  
  # estimate Highest Posterior Interval
  hh = Hpi(df, binned = TRUE)*1
  
  # make 2D dimensional smoothing
  dd = kde(df, H = hh)
  
  # create contour lines for polygons
  cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                    z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                  length = 20))
  
  # set color
  col1 = alpha(brewer.pal(3, col.pals[4]), 0.4)
  
  # plot three stacked polygons for each level of the HPI
  polygon(cc[[16]]$x, cc[[16]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
  polygon(cc[[17]]$x, cc[[17]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
  polygon(cc[[18]]$x, cc[[18]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
}



################
# Each tribe
################
#color palettes in order of polygons
col.pals = c('Reds', 'Greens', 'Oranges', 'Greys', 'Purples','BuPu','Blues')
################
# Ex: Akodontini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 374 é Akodontini
mrca(tree.sigmodontinae)["Akodon_azarae", "Kunsia_tomentosus"] # 374
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 114
# plot map
plot(shape.neot)
# Loop for plotting the one polygon

# subset of rase results
rase_results_subset<-rase_results[2001:12000,]

# make `x` & `y` object
df = data.frame(rase_results_subset[,114], rase_results_subset[,114+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
mean(rase_results_subset[,114])
mean(rase_results_subset[,114+(ncol(rase_results_subset)/2)-1])
points(-58.76,-21.12)

#bifurcating 1 #N 375
df = data.frame(rase_results_subset[,115], rase_results_subset[,115+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
#again
df = data.frame(rase_results_subset[,126], rase_results_subset[,115+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
df = data.frame(rase_results_subset[,116], rase_results_subset[,116+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
#bifurcating 2 #N 429
df = data.frame(rase_results_subset[,169], rase_results_subset[,169+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
#again
df = data.frame(rase_results_subset[,170], rase_results_subset[,170+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
df = data.frame(rase_results_subset[,176], rase_results_subset[,176+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = alpha(brewer.pal(3, col.pals[1]), 0.4)

# plot three stacked polygons for each level of the HPI
# 20 classes, então 5% pra cada. A 19 é 95%, a 18 é 90%, assim por diante (15=75%; 10=50%)
# 89 max, então 70 min - 95% = 88, 90%=87, 80%=85, 75%=84, 70%=83, 50%=79
polygon(cc[[87]]$x, cc[[87]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[85]]$x, cc[[85]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[83]]$x, cc[[83]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[87]]$y)
range(cc[[87]]$y)
mean(cc[[87]]$x)
range(cc[[87]]$x)

#bifurcating 1 #N 375
polygon(cc[[56]]$x, cc[[56]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[54]]$x, cc[[54]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[52]]$x, cc[[52]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
#again
polygon(cc[[53]]$x, cc[[53]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[51]]$x, cc[[51]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[49]]$x, cc[[49]]$y, col = col1[3], border = 'grey60', lwd = 0.3)

polygon(cc[[136]]$x, cc[[136]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[134]]$x, cc[[134]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[132]]$x, cc[[132]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
#bifurcating 2 #N 429
polygon(cc[[85]]$x, cc[[85]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[83]]$x, cc[[83]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[81]]$x, cc[[81]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
#again
polygon(cc[[82]]$x, cc[[82]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[80]]$x, cc[[80]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[78]]$x, cc[[78]]$y, col = col1[3], border = 'grey60', lwd = 0.3)

polygon(cc[[90]]$x, cc[[90]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[88]]$x, cc[[88]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[86]]$x, cc[[86]]$y, col = col1[3], border = 'grey60', lwd = 0.3)



# Ex: Oryzomyini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 438 é Oryzomyini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 178
# plot map
plot(shape.neot)
# Loop for plotting the one polygon

# subset of rase results
rase_results_subset<-rase_results[2001:12000,]

# make `x` & `y` object
df = data.frame(rase_results_subset[,178], rase_results_subset[,178+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
mean(rase_results_subset[,178])

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = alpha(brewer.pal(3, col.pals[2]), 0.4)

# plot three stacked polygons for each level of the HPI
# 63 max, 44 min - 62=95%, 61=90%, 80%=59, 70%=57
polygon(cc[[61]]$x, cc[[61]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[59]]$x, cc[[59]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[57]]$x, cc[[57]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[61]]$y)
range(cc[[61]]$y)
mean(cc[[61]]$x)
range(cc[[61]]$x)

# Ex: Phyllotini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 279 é Phyllotini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 19
# plot map
plot(shape.neot)
# Loop for plotting the one polygon

# make `x` & `y` object
df = data.frame(rase_results_subset[,19], rase_results_subset[,19+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = alpha(brewer.pal(3, col.pals[3]), 0.4)

# plot three stacked polygons for each level of the HPI
# 63 max - 61=90%, 59=80%, 57=70%
polygon(cc[[61]]$x, cc[[61]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[59]]$x, cc[[59]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[57]]$x, cc[[57]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[61]]$y)
range(cc[[61]]$y)
mean(cc[[61]]$x)
range(cc[[61]]$x)

# Ex: Thomasomyini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 347 é Thomasomyini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 87
# plot map
plot(shape.neot)
# Loop for plotting the one polygon

# make `x` & `y` object
df = data.frame(rase_results_subset[,87], rase_results_subset[,87+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = alpha(brewer.pal(3, col.pals[4]), 0.4)

# plot three stacked polygons for each level of the HPI
# 80 max, 78=90%, 76=80%, 74=70%
polygon(cc[[78]]$x, cc[[78]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[76]]$x, cc[[76]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[74]]$x, cc[[74]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[78]]$y)
range(cc[[78]]$y)
mean(cc[[78]]$x)
range(cc[[78]]$x)

# Ex: Abrotrichini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 326 é Abrotrichini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 66
# plot map
plot(shape.neot.beh)
# Loop for plotting the one polygon

# make `x` & `y` object
df = data.frame(rase_results_subset[,66], rase_results_subset[,66+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
mean(rase_results_subset[,66])

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = addalpha(brewer.pal(3, col.pals[5]), 0.4)

# plot three stacked polygons for each level of the HPI
# 66 max - 64=90%, 62=80%, 60=70%
polygon(cc[[64]]$x, cc[[64]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[62]]$x, cc[[62]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[60]]$x, cc[[60]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[64]]$y)
range(cc[[64]]$y)
mean(cc[[64]]$x)
range(cc[[64]]$x)



# Ex: Sigmodontini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 263 é Sigmodontini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 3
# plot map
plot(shape.neot.beh)
# Loop for plotting the one polygon

# make `x` & `y` object
df = data.frame(rase_results_subset[,3], rase_results_subset[,3+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
mean(rase_results_subset[,3])
mean(rase_results_subset[,3+(ncol(rase_results_subset)/2)-1])

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = addalpha(brewer.pal(3, col.pals[6]), 0.4)

# plot three stacked polygons for each level of the HPI
cc
# 97 max, 95=90%, 93=80%, 91=70%
polygon(cc[[95]]$x, cc[[95]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[93]]$x, cc[[93]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[91]]$x, cc[[91]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[95]]$y)
range(cc[[95]]$y)
mean(cc[[95]]$x)
range(cc[[95]]$x)



# Ex: Ichtyomyini
plot(tree.sigmodontinae,cex=0.3)
nodelabels(cex=0.3) # nó 273 é Sigmodontini
as.data.frame(dimnames(rase_results)[2]) # corresponde à recon. N 13
# plot map
plot(shape.neot.beh)
# Loop for plotting the one polygon

# make `x` & `y` object
df = data.frame(rase_results_subset[,13], rase_results_subset[,13+(ncol(rase_results_subset)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
mean(rase_results_subset[,13])
mean(rase_results_subset[,13+(ncol(rase_results_subset)/2)-1])

# estimate Highest Posterior Interval
hh = Hpi(df, binned = TRUE)*1

# make 2D dimensional smoothing
dd = kde(df, H = hh)

# create contour lines for polygons
cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                  z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                length = 20))

# set color
col1 = addalpha(brewer.pal(3, col.pals[7]), 0.4)

# plot three stacked polygons for each level of the HPI
cc
# 85 max, 83=90%, 81=80%, 79=70%
polygon(cc[[83]]$x, cc[[83]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
polygon(cc[[81]]$x, cc[[81]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
polygon(cc[[79]]$x, cc[[79]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
mean(cc[[83]]$y)
range(cc[[83]]$y)
mean(cc[[83]]$x)
range(cc[[83]]$x)  




# Loop for plotting all polygon for nodes
for (i in 1:(ncol(rase_results)/2)) {
  
  # make `x` & `y` object
  df = data.frame(rase_results[,i], rase_results[,i+(ncol(rase_results)/2)-1]) # colunas vezes 2 - sigmax e y (2 colunas (-1)), pra achar a posição do y
  
  # estimate Highest Posterior Interval
  hh = Hpi(df, binned = TRUE)*1
  
  # make 2D dimensional smoothing
  dd = kde(df, H = hh)
  
  # create contour lines for polygons
  cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]], 
                    z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                  length = 20))
  
  # set color
  #col1 = addalpha(brewer.pal(3, col.pals[4]), 0.4)
  
  # plot three stacked polygons for each level of the HPI
  polygon(cc[[16]]$x, cc[[16]]$y, col = col1[1], border = 'grey60', lwd = 0.3) 
  polygon(cc[[17]]$x, cc[[17]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
  polygon(cc[[18]]$x, cc[[18]]$y, col = col1[3], border = 'grey60', lwd = 0.3)
}


