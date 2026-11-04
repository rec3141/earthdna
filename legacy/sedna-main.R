#run sedna-cluster-seed.sh to generate the 97% clusters

#/opt/local/bin/R
#.libPaths("/Users/cryomics/Library/MRO/3.3/library") #if using MRO

library(plyr)
library(data.table)
library(smacof)
library(rgl)
library(phangorn)
library(sphereplot)
library(phytools)
library(spatstat)
library(parallel)


setwd("/work/cryomics/reference_dbs/SILVA/sedna")
silva.acc_taxid <- read.table("./../tax_slv_ssu_123.acc_taxid",row.names=1, stringsAsFactors=FALSE)
silva.acc_taxid[,2] <- sapply(rownames(silva.acc_taxid),function(x) strsplit(x,"[.]")[[1]][1])
colnames(silva.acc_taxid) <- c("taxid","acc")

silva.map <- read.table("./../tax_slv_ssu_123.map",stringsAsFactors=FALSE,sep="\t",head=T)
silva.map <- silva.map[,-(3:4)]
colnames(silva.map) <- c("taxid","taxa")

silva.tat <- join(silva.acc_taxid,silva.map)

silva.tax <- read.table("./../silva.seed_v123.full.tax",stringsAsFactors=FALSE,sep=";",head=F,fill=T)
silva.tax <- cbind(matrix(unlist(strsplit(silva.tax$V1,split="\t")),ncol=2,byrow=T),silva.tax)
rownames(silva.tax) <- silva.tax[,1]
silva.tax <- silva.tax[,-c(1,3,9)]
colnames(silva.tax) <- taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus") #infrakingdom missing... i dunno

dist.rows <- as.numeric(strsplit(system("wc -l ./../silva.seed_v123.filter.square.dist",intern=T)," +")[[1]][2])
dist.read <- fread("./../silva.seed_v123.filter.square.dist",skip=1,stringsAsFactors=F,nrows=dist.rows)
dist.labels <- dist.read$V1

#need to trim read names and final whitespace
dist.in <- dist.read[,-c(1,ncol(dist.read)),with=F]
dist.otus.in <- read.table("./../silva.seed_v123.filter.square.an.0.03.rep.names",sep="\t",stringsAsFactors=F)
dist.select <- dist.labels %in% dist.otus.in$V1

#another alternative: interpolate parent nodes a la Indiana, then heavily weight those nodes in MDS

#have a tree alongside, clicking node on tree highlights group on map

#for testing
#dist.select <- dist.labels %in% dist.otus.in$V1[seq(1,6900,24)]
#dist.select <- as.logical(dist.labels %in% dist.otus.in$V1[seq(1,6900,24)+1] + dist.select)
dist.mat <- as.matrix(dist.in[dist.select,dist.select,with=F])
rownames(dist.mat) <- dist.labels[dist.select]
colnames(dist.mat) <- dist.labels[dist.select]

dist.weight <- as.dist(dist.mat^.5)*1

#start at 19:24
date()
dist.sphere5 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=T,modulus=1,init="random")
#OP3 spread way out in halo, don't like shape of bacteria
date()
dist.sphere6 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=F,modulus=1,init="random")
#messy
date()
dist.sphere1 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=T,modulus=1)
date()
dist.sphere2 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=F,modulus=1)
date()
dist.sphere3 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=T,modulus=2)
date()
dist.sphere4 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500,penalty=2e3,verbose=T,eps=1,weightmat=dist.weight,type="ratio",relax=F,modulus=2)
date()
dist.sphere7 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=5000,penalty=2e3,verbose=T,eps=0.1,weightmat=dist.weight,type="ratio")
date()
dist.sphere8 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=5000,penalty=2e3,verbose=T,eps=0.01,weightmat=dist.weight,type="ratio")
date()
dist.sphere9 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=50000,penalty=2e3,verbose=T,eps=0.1,weightmat=dist.weight,type="ratio")
date()
dist.sphere10 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=50000,penalty=2e3,verbose=T,eps=0.01,weightmat=dist.weight,type="ratio")
date()
dist.sphere11 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=50000,penalty=2e3,verbose=T,eps=0.001,weightmat=dist.weight,type="ratio")
date()
dist.sphere12 <- smacofSphere(dist.mat,"dual", ndim=3, itmax=500000,penalty=2e3,verbose=T,eps=0.0001,weightmat=dist.weight,type="ratio")
date()


#save(file="dist.sphere.0.03.1e4_1.save.R",dist.sphere) #2160 iterations stress 1190679.35159630
#save(file="dist.sphere_sq.0.03.1e4_1.save.R",dist.sphere_sq) #Iteration:  2409  Stress (not normalized):  366675.32561225
#load("dist.sphere.0.03.1e4_1.save.R")

dist.sphere <- dist.sphere_wt
dist.sphere <- dist.sphere_grampos
dist.sphere <- dist.sphere_gammas
dist.sphere <- dist.sphere_gammas_wt12

dist.sphere <- dist.sphere_all_wt17

dist.sphere <- dist.sphere1

### now plot it
silva97.conf <- dist.sphere$conf
rownames(dist.sphere$conf) <- rownames(dist.mat)

#normalize to earth radius
radius.earth <- 6378137
silva97.sph <- car2sph(silva97.conf)
silva97.sph[,3] <- radius.earth
rownames(silva97.sph) <- rownames(dist.mat)
silva97.xyz <- sph2car(silva97.sph)
open3d(windowRect=c(0,0,500,500))
bg3d(sphere=F, back = "filled", color = "lightgrey")
spheres3d(0,0,0, radius = 0.97*radius.earth, color = "lightblue", alpha = 1, back = "fill", front = "fill", lit=T, specular="black",  ax.grid = TRUE, sphere.rgl = TRUE)
#rgl.sphtext(silva97.sph, text=silva.tax[rownames(silva97.sph),"order"],cex=0.5)
#rgl.sphtext(silva97.sph, text=rownames(silva97.sph))
rgl.sphtext(silva97.sph, text=".")




write.table(cbind(silva97.sph,silva.tax[rownames(silva97.sph),]),file="dist.sphere_sq.sph.csv",sep="\t",quote=F,col.names=1)

## 100 reps ratio distx10,000 pen1e7 wt^(1/3): orders not well defined
## 100 reps ratio distx10,000 pen1e7 wt^(1/2): phyla tightly clustered, no good



##### Now rotate so that we have nicely distributed continents

#functions to convert cartesian to radians and degrees

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

#where do we want the new north pole? long, lat, radius
georef.deg <- t(as.matrix(c(90, 15, radius.earth))) #empirically determined, will need to change for other datasets
georef.x <- deg2rad(12) #empirically determined
georef.car <- sph2car(georef.deg,deg=T)
georef.rad <- car2sph(georef.car,deg=F)

#rotate in z dimension (longitude) to the prime meridian
silva97.rotz.xyz <- asEuclidean(rotate3d(asHomogeneous(silva97.xyz), georef.rad[1], 0, 0, 1))
silva97.rotz.sph <- car2sph(silva97.rotz.xyz, deg=T)

#rotate in y dimension to ...
silva97.roty.xyz <- asEuclidean(rotate3d(asHomogeneous(silva97.rotz.xyz), pi/2 - georef.rad[2], 0, 1, 0))
silva97.roty.sph <- car2sph(silva97.roty.xyz, deg=T)

#rotate in x dimension to ...
silva97.rotx.xyz <- asEuclidean(rotate3d(asHomogeneous(silva97.roty.xyz), georef.x, 0, 0, 1))
silva97.rotx.sph <- car2sph(silva97.rotx.xyz, deg=T)
open3d(windowRect=c(0,0,800,800))
bg3d(sphere=F, back = "filled", color = "lightgrey")
spheres3d(0,0,0, radius = 0.97*radius.earth, color = "lightblue", alpha = 1, back = "fill", front = "fill", lit=T, specular="black",  ax.grid = TRUE, sphere.rgl = TRUE)
rgl.sphgrid(radius = 1.02*radius.earth, col.long='red', col.lat='blue', deggap = 15, longtype = "D", add = TRUE, radaxis=FALSE)
rgl.sphtext(silva97.rotx.sph,text=".",col="blue")

silva97.map.sph <- silva97.rotx.sph
silva97.map.xyz <- sph2car(silva97.rotx.sph,deg=T)
silva97.map.rad <- car2sph(silva97.map.xyz, deg=F)

plot(silva97.map.sph[,1:2])

### rearranging each sequence to north pole, then finding latitude enclosing 128 closest sequences

high.lat.save <- list()
#georef is row to place at north pole
for (georef in 1:nrow(silva97.map.xyz)) {
	print(georef)
	seq.rad <- silva97.map.rad 
	
	#rotate in z dimension (longitude) to the prime meridian
	car.rot.lon <- asEuclidean(rotate3d(asHomogeneous(silva97.map.xyz), seq.rad[georef,1], 0, 0, 1))
	sph.rot.lon <- car2sph(car.rot.lon, deg=T)

	seq.rad <- car2sph(car.rot.lon, deg=F)
	colnames(car.rot.lon) <- c("x","y","z")
	car.rot.lat <- asEuclidean(rotate3d(asHomogeneous(car.rot.lon), pi/2 - seq.rad[georef,2], 0, 1, 0))
	sph.rot.lat <- car2sph(car.rot.lat, deg=T)

	# row 'georef' is now at the north pole

	#get the highest latitute 128 points + georef
	high.lat.save[[georef]] <- rownames(sph.rot.lat[order(sph.rot.lat[,2],decreasing=T),][2:129,])

}

high.lat.write <- t(as.data.frame(high.lat.save))
rownames(high.lat.write) <- rownames(silva97.map.xyz)

write.table(file="high.lat.out.csv",high.lat.write,sep="\t",quote=F,row.names=T,col.names=F)
write.table(file="sampleMappingFile.silva.txt",silva97.map.xyz,quote=F,sep="\t",row.names=T,col.names=F)


##### DONE HERE


### we want to interpolate NR99, then use those distances + silva97 to interpolate user seqs
#!/usr/bin/env bash

bash sedna-interpolate-nr.sh





#!/usr/bin/env bash

#FOR USER

#classify against SEED? otherwise will have to classify NR hit against SEED to map it

#awk '{ if ($16>97) a[$3]++ } END {for(i in a) if(a[i]>1) print i,a[i]}' $1 | awk -F'[. ]' '{print $1,$3}' > tmp.me.weights
#awk '{ if ($2>100) print $3}' ambon.trim.contigs.unique.align.report
 
#do more here for <97% similar

# do the interpolation

bash sedna-interpolate-nr.sh &
#or
bash mdsearth-run-mothur.sh &
#or
bash sedna-run-mothur-dada.sh &


#use NLS interpolation?





#!/opt/local/bin/R
library(sp)
#library(parallel)
library(doMC)
library(foreach)
#... others above
library(alphahull)
library(plyr)
library(Matrix)
library(rgdal)
library(rgeos)
library(ape)
library(phangorn)
# installation on MacOSX is stupid, following instructions here: http://tlocoh.r-forge.r-project.org/mac_rgeos_rgdal.html
# 1. download source e.g. http://cran.at.r-project.org/src/contrib/rgeos_0.3-19.tar.gz
# 2. do [sudo?] /opt/local/bin/R CMD INSTALL rgeos_0.3-19.tar.gz --configure-args="--with-geos-config=/Library/Frameworks/GEOS.framework/unix/bin/geos-config"

library(raster)
library(sphereplot)

#load all interpolated files
indir <- "nr"
infile <- "nr_all.mds" #pre-catted: much faster than looping
mds.save <- read.table(paste0(indir,"/",infile),stringsAsFactors=F)

mds.in <- mds.save  #make a copy to work with
#count.mds.in <- count(cbind(rownames(mds.save),mds.save))  # counts number of times each reference sequence was seen
mds.in <- unique(mds.in) #do unique before dropping names, otherwise lose 100,000 seqs
rownames(mds.in) <- mds.in[,1]
mds.in <- mds.in[,-1]
count.mds.in <- count(mds.in)
count.mds.all <- join(mds.in,count.mds.in)
rownames(count.mds.all) <- rownames(mds.in)

radius.earth <- 6378137
pick.mds.sample <- mds.in[,4]>0
pick.mds.silva97 <- mds.in[,4]==0
mds.in <- mds.in[,-4]
mds.sph <- car2sph(mds.in)
mds.sph[,3] <- radius.earth
rownames(mds.sph) <- rownames(mds.in)

write.table(file="nr_all.mds.sph",mds.sph,col.names=T,quote=F)

writeOGR(silva97.map.sph,"silva97","silva97",driver="ESRI Shapefile",overwrite=T)
writeOGR(mds.sph,"mds_sph","mds_sph",driver="ESRI Shapefile",overwrite=T)




### start mapping NR

# NOT AUTOMATED
# in ARB, reroot tree and add acc_name NDS, export tree
# in ARB, add a new field called acc_name, then export tree using that via NDS
#   :*=*(acc).*(name)
# called it SSURefNR99_1200_slv_123_name_acc.ntree
# make sure to set NDS width higher than acc_name width so they aren't cropped
# nodetype: nds; save branch lengths: yes; save bootstrap values: no; save group names: yes; name quoting: no; replace problem chars: yes;
# remove header inserted by ARB
# remove newlines: tr -d $'\n' < SSURefNR99_1200_slv_123_name_acc.ntree > SSURefNR99_1200_slv_123_name_acc.tree

#cleaned up: tr -d $'\n' < SSURefNR99_1200_slv_123_clean.ntree > SSURefNR99_1200_slv_123_clean.tree
#edited tree to remove long branches to clean up labels
#1_ 10% 0.001
#2_ 50 0.005
#3_ 90 0.01
#4_ 50 0.1
#5_ remove red
#6_ replace blue and seed97

#it's a big tree, this takes a minute; read.newick didn't finish
#tree.arb <- read.tree(file="SSURefNR99_1200_slv_123_name_acc.tree")

tree.arb.clean <- read.tree(file="SSURefNR99_1200_slv_123_clean.tree")
tree.arb <- tree.arb.clean

nc <- 14 #number of cores
registerDoMC(cores=nc)

source("sort.edges.r")
source("get.polygons.r")

### first do nodes with labels

#using nodes in tree with labels
node.list <- tree.arb$Nnode + 1 + which(tree.arb$node.label != "")
#makes a list of all the polygons
map.labeled <- foreach(tree.node=node.list,.inorder=F,.errorhandling='remove') %dopar% get.polygons(tree.node, tree.arb, mds.sph)
#turns it into a SpatialPolygons object
map.sp.labeled <- SpatialPolygons(map.labeled)
map.sp.taxa <- as.numeric(unlist(lapply(map.sp.labeled@polygons,function(x) x@ID)))
map.sp.edges <- tree.arb$edge.length[map.sp.taxa]


### now output ESRI shapefiles
map.spdf.row.names <- unlist(lapply(map.sp.labeled@polygons,function(x) x@ID))
map.spdf.elevation <- rep(1,length(map.sp.labeled))
map.spdf.taxonomy <- tree.arb$node.label[as.numeric(map.spdf.row.names) - tree.arb$Nnode - 1]
map.spdf.labeled <- SpatialPolygonsDataFrame(map.sp.labeled,data=data.frame("elevation"=map.spdf.elevation,"taxonomy"=map.spdf.taxonomy,row.names=map.spdf.row.names))
writeOGR(map.spdf.labeled,"nr_labels_shapes","nr_labels_shapes",driver="ESRI Shapefile",overwrite=T)

map.sl.labeled <- as(map.spdf.labeled, 'SpatialLinesDataFrame')
writeOGR(map.sl.labeled,"nr_labels_lines","nr_labels_lines",driver="ESRI Shapefile",overwrite=T)

### now rasterize

map.nodes <- as.numeric(map.spdf.row.names)
map.areas <- as.numeric(unlist(lapply(map.labeled,function(x) x@area)))

#weights as number of descendants

#weights.descent <- unlist(mclapply(map.nodes, function(tree.node) length(unlist(Descendants(tree.arb, tree.node, type="tips"))), mc.cores=14))
#weights.raster <- weights.descent

#weights as sqrt(number of descendants)
#weights.raster <- sqrt(weights.descent)

#weights as number of descendants/area
weights.raster <- weights.descent/map.areas
weights.raster[is.infinite(weights.raster)] <- 0

#weights as log(number of descendands/area)
#weights.raster <- log(weights.descent/map.areas)
#weights.raster[is.infinite(weights.raster)] <- 0

#weights as positive log(number of descendands/area)
weights.raster <- log(weights.descent/map.areas)+abs(min(log(weights.descent/map.areas)))
weights.raster[is.infinite(weights.raster)] <- 0

#weights as 1
#weights.raster <- rep(1,length(map.sp.labeled))
#weights as ege length
#weights.raster <- map.sp.edges

blank.raster <- raster(ncol=720, nrow=360,xmn=-180,xmx=180,ymn=-90,ymx=90)
x<-extent(blank.raster)[1]:extent(blank.raster)[2]
numchunks <- nc*6
chunks <- split(x, sort(x%%numchunks))

raster.comb <- function(a,b) merge(a,b)

map.strips <- foreach(strip=1:length(chunks),.combine='raster.comb',.inorder=T,.errorhandling='remove') %dopar% {
	xmin <- min(chunks[[strip]])
	if(strip == length(chunks)) { xmax <- max(chunks[[strip]]) } else {xmax <- max(chunks[[strip]]) + 1 }
	ymin <- extent(blank.raster)[3]
	ymax <- extent(blank.raster)[4]
	
	rasterize(map.sp.labeled, crop(blank.raster, extent(xmin, xmax, ymin, ymax) ), field=weights.raster, fun='sum',progress='text')
} 

plot(map.strips)
writeRaster(map.strips, file="map.labeled3.tif",format="GTiff",overwrite=T)


### then do it for everything

#using nodes in tree with labels
node.list <- seq(length(tree.arb$tip.label)+1,length(tree.arb$edge.length)+1)
#makes a list of all the polygons
map.nrclean <- foreach(tree.node=node.list,.inorder=F,.errorhandling='remove') %dopar% get.polygons(tree.node, tree.arb, mds.sph)







#everything

#using all nodes in tree, looking like 12-16 hours
node.list <- seq(length(tree.arb$tip.label)+1,length(tree.arb$edge.length)+1)
#makes a list of all the polygons
map.polygons <- foreach(tree.node=node.list,.inorder=F,.errorhandling='remove') %dopar% get.polygons(tree.node, tree.arb, mds.sph)
#turns it into a SpatialPolygons object
map.sp <- SpatialPolygons(map.polygons)

#makes blank canvas
blank.raster <- raster(ncol=7200, nrow=3600,xmn=-180,xmx=180,ymn=-90,ymx=90)
weights.raster <- rep(1,length(map.plot)) #can change this to, e.g. edge.length

map.raster <- rasterize(map.sp, blank.raster, field=weights.raster,fun='sum',progress='text')
writeRaster(map.raster, file="map.raster.tif",format="GTiff",overwrite=T)




blank.raster <- raster(ncol=720, nrow=360,xmn=-180,xmx=180,ymn=-90,ymx=90)
map.raster.labeled <- rasterize(map.sp.labeled, blank.raster, field=rep(1,length(map.sp.labeled)),fun='sum',progress='text')
writeRaster(map.raster, file="map.labeled.tif",format="GTiff",overwrite=T)









#IMPORT EACH SAMPLE

### 20 samples in FAST mode, 9:35-10:19 = 45 minutes
### 812 samples in SLOW mode, 6 hours
### calculating distance matrix from tree in ARB: used 30Gb RAM, still 0% complete after 30minutes

#plot each sample and abundance

#import count table: rows are unique seqs, cols are samples
#for each sample we can plot which sequences and how many
count.table <- read.table("./NOAA-testing/input.trim.contigs.pcr.subsample.count_table",head=T,row.names=1)
#dada gives this (check orientation)

count.table.names <- unlist(lapply(rownames(count.table),function(x) substr(x,1,nchar(x)-1)))
rownames(count.table) <- count.table.names

###need to add uniques or did mothur do it for us? mothur did it for us -- how nice
#count.join <- merge(as.data.frame(mds.sph),count.table)
colors <- rainbow(ncol(count.table))

for (pick.sample in c(2,20,220,320,400,405,425,450)) {
	print(colnames(count.table)[pick.sample])
	pick.mds.single <- count.table.names[count.table[,pick.sample]>0]
	mds.single <- mds.save[pick.mds.single,]
	mds.single.sph <- car2sph(mds.single[,1:3],deg=T)
	mds.single[,3] <- radius.earth*1.1
#	rgl.sphtext(mds.single.sph,text='.',cex=count.table[pick.mds.single,5], color=colors(pick.sample))
	rgl.sphtext(mds.single.sph,text='.',cex=3, color=colors(pick.sample))
}

###### CODE FOR TURNING RELATIVE ABUNDANCE INTO ALTITUDE
### BASIC IDEA: 
### each node in the tree defines a convex hull polygon on the surface of the sphere
### step 1: define polygons
### step 2: for each sequence, find which polygons it sits on
### step 3: for each polygon, add an altitude proportional to the number of reads in the sequence
### 	will need to try various schemes for allocating altitude so that it is visually appealing
###		method 1: height added to each polygon = number of reads in sequence -- probably ugly
###	*	method 2: height added to each polygon = sqrt(number of reads in sequence)
###		method 3: height added to each polygon = number or sqrt(number) of reads in all sequences within polygon
###		method 4: height added to each polygon = number or sqrt(number) of reads scaled by area
###	*	method 5: height added to each polygon = number or sqrt(number) of reads scaled by branch length on tree

#R
library(smacof)
library(rgl)
library(phangorn)
library(sphereplot)
library(phytools)
library(spatstat)

#system("ln -s /work/cryomics/reference_dbs/SILVA/SSURefNR99_1200_slv_123.ntree")
#silva.tree <- read.tree("SSURefNR99_1200_slv_123.ntree")
# v123.1 has different names from v123
#use SEED 97% as basic platform --update to use NR99 at some point, now using mds.in intermediate
#seq.arb <- cbind(as.character(sapply(rownames(silva97.xyz),function(x) strsplit(x,'[.]')[[1]][1])),as.character(sapply(rownames(silva97.xyz),function(x) strsplit(x,'[.]')[[1]][2])))
#seq.accs <- unique(sort(as.character(sapply(rownames(silva97.xyz),function(x) strsplit(x,'[.]')[[1]][1]))))


seq.names <-  unique(sort(as.character(sapply(rownames(mds.in),function(x) strsplit(x,'[.]')[[1]][2]))))
tree.arb <- stree(length(seq.names), type = "star", tip.label = seq.names)
write.tree(tree.arb,file="arb.out.tree")

#NOT AUTOMATED
#eventually want to use entire NR99 tree here...
#open in arb SSURef_NR99_123_SILVA_12_07_15_opt.arb, export reference tree
#open in arb SSURef_123_SILVA_19_07_15_opt.arb, import reference tree
#import arb.out.tree
#mark species in arb.out.tree
#select reference tree
#remove species not marked
#reroot to (B,(A,E))
#export pruned reference tree as SSURefNR99_1200_slv_123_trim.tree
#remove header inserted by ARB
#remove newlines: tr -d $'\n' < SSURefNR99_1200_slv_123.reroot.tree > SSURefNR99_1200_slv_123.reroot.nr.tree
#seq.tree.arb <- read.newick(file="SSURefNR99_1200_slv_123_trim_nr.tree")
seq.tree.arb <- read.newick(file="SSURefNR99_1200_slv_123.reroot.nr.tree")
#if we reorient world to put prime meridian between bacteria and archaea then we can do area calcs on 2d and not worry about edges (?)
#not necessarily because shortest line in stereo is not same as shortest line in ortho (??) but let's do it, it's faster and easier

#adding here to use mds.in instead of silva97.xyz, this might break something to do with the tree

#### moved nice rotation up up up
#need this anymore?
mds.xyz <- as.matrix(mds.in)
mds.lle <- car2sph(mds.xyz,deg=T)
mds.lle[,3] <- radius.earth
mds.xyz <- sph2car(mds.lle,deg=T)
colnames(mds.xyz) <- c("x","y","z")

#mercator plot...
#silva97.map.sphy <- cbind(silva97.map.sph[,1],silva97.map.sph[,2]) #cut out conversion to 0->360

#now let's do 2d convex hull on map
#or... just do it by tree, which is going to look better?
#problem now where convex hull contains many points that are not in the subtree
#ahh but that's the point, we want the interpolated/user sequences to get heights
#but if we just do convex hull they'll be double counted
#may need to use pplacer to put them on the tree... ?

library(sp)
library(raster)
library(rgdal)
library(Matrix) #for sparse matrices
library(foreach) #for parallel
library(doMC) #multicore
library(alphahull)

# install.packages('rgdal', type = "source", configure.args=c(
#     '--with-gdal-config=/Library/Frameworks/GDAL.framework/Programs/gdal-config',
#     '--with-proj-include=/Library/Frameworks/PROJ.framework/Headers',
#     '--with-proj-lib=/Library/Frameworks/PROJ.framework/unix/lib'))
     
#using silva97.xyz here because we want only sequences in the tree
sph.names <- as.character(sapply(rownames(silva97.map.sph),function(x) strsplit(x,'[.]')[[1]][2]))
#seq.names <- as.character(sapply(rownames(silva97.xyz),function(x) strsplit(x,'[.]')[[1]][2]))
#seq.deg <- car2sph(silva97.xyz, deg=T)

map.alt <- rep(0,nrow(silva97.map.sph))
map.poly <- list()
map.size.lon <- 7201
map.size.lat <- ((map.size.lon-1)/2)+1
grid.lat <- matrix(data=1,nrow=map.size.lat,ncol=map.size.lon)*seq(from=90,to=-90,length.out=map.size.lat)
grid.lon <- t(matrix(data=1,nrow=map.size.lat,ncol=map.size.lon)*seq(from=-180,to=180,length.out=map.size.lat))

# if the tree doesn't change...
# we should only need to calculate pip's once and then apply that to all of the samples to quickly make surfaces for each
# since each tree node is independent we should be able to split this up into subtrees and run them simultaneously

#source("treeSlice.R") #cuts tree, returns subtrees
# we're going to calculate each pip and store it in a sparse matrix
# row=lat, col=lon, z= [1,0] if that point is in the polygon
# we'll repeat this for each node of the tree

#I think I should make a stacked vector array rather than the raster
#then do the rastering within GDAL or MapBox or ..

registerDoMC(cores=14)

#function to get pip
source("get.pip.grid.r")
source("sort.edges.r")
source("get.pip.seqs")
source("sparse.sum.r")





# this makes the grid for each node of the tree
map.layers <- list()
map.layers <- foreach(tree.node=(length(tree.arb$tip.label)+1):nrow(seq.tree.arb$edge)) %dopar% get.pip.grid(tree.node, seq.tree.arb, sph.names, silva97.map.sph, grid.lon, grid.lat, map.size.lat, map.size.lon)

#for(tree.node in (length(tree.arb$tip.label)+1):nrow(seq.tree.arb$edge)) get.pip.grid(tree.node, seq.tree.arb, sph.names, silva97.map.sph, grid.lon, grid.lat, map.size.lat, map.size.lon)

map.grid <- sparse.sum(map.layers)
map.raster <- raster(as.matrix(map.grid),xmn=-180,xmx=180,ymn=-90,ymx=90)
writeRaster(map.raster, "map.export.geo", format="GTiff",overwrite=T)

map.out <- cbind(silva97.map.sph[,1:2],map.alt) #map.alt doesn't exist any more... now we want count.table ? or rowSums(something)
map.out <- cbind(rownames(map.out),map.out)
colnames(map.out) <- c("name","longitude","latitude","elevation")
map.export <- cbind(map.out,silva.tax[rownames(map.out),])
map.export <- map.export[,-ncol(map.export)]
i <- sapply(map.export, is.factor)
map.export[i] <- lapply(map.export[i], as.character)

write.table(file=paste0("unknown.export.csv"),map.export[which(is.na(map.export$domain)),],sep="\t",quote=F,row.names=F)

write.table(file=paste0("known.export.csv"),map.export[which(!is.na(map.export$domain)),],sep="\t",quote=F,row.names=F)

seq.pips <- list()
seq.pips <- foreach(tree.node=(length(tree.arb$tip.label)+1):nrow(seq.tree.arb$edge)) %dopar% get.pip.seqs(tree.node, seq.tree.arb, sph.names, silva97.map.sph, grid.lon, grid.lat, map.size.lat, map.size.lon)
mat.pips <- do.call(cbind, seq.pips) #this gives rows=seqs, cols=nodes

#ALMOST THERE

#ok first try let's add the total number of sequences in a node to that layer
#first pick a sample
#second find which rows are in that sample
#third for each layer, find which of those rows are present
#fourth sum the rows and multiply rows

sph.rot.names <- rownames(silva97.map.sph)

layer.seqs <- Matrix(0,ncol=ncol(count.table),nrow=ncol(mat.pips)) #rows=nodes, cols=samples
for(i in 1:ncol(count.table)) {
	print(i)
	layer.sums <- foreach(j=1:ncol(mat.pips),.combine=c) %dopar% sum(count.table[sph.rot.names[mat.pips[,j]>0],i],na.rm=T)
	layer.seqs[,i] <- layer.sums
}

#ok now layer.seqs hold the total number of sequences in a node in each sample
#rows=nodes, cols=samples
#we want to multiply each layer by that number and sum them up into a topography

#more parallel version
sample.rasters <- foreach(i=1:ncol(layer.seqs)) %dopar% {
	sample.layers <- list()
	for (j in 1:nrow(layer.seqs)) {
	#this is where we can change the weighting
		sample.layers[[j]] <- map.layers[[j]] * sqrt(layer.seqs[j,i]) * seq.tree.arb$edge.length[j+1+seq.tree.arb$Nnode] #length(tree.arb$tip.label)+1
#		sample.layers[[j]] <- map.layers[[j]] * sqrt(layer.seqs[j,i])		
#		sample.layers[[j]] <- map.layers[[j]] * layer.seqs[j,i]		
	}
	map.samples <- sparse.sum(sample.layers)
	raster(as.matrix(map.samples),xmn=-180,xmx=180,ymn=-90,ymx=90)
}

#print them out
for(i in 1:length(sample.rasters)) {

	writeRaster(sample.rasters[[i]], paste0("sample.raster.sqrt_edge.",i,".tif"), format="GTiff",overwrite=T)
	
	sample.select.x <- rownames(silva97.map.sph) %in% count.table.names[count.table[,i]>0]
	#these aren't reciprocal, not clear why at this point but we'll work around it
	#there are a bunch of sequences in count.table but not silva97.map.sph ... why?
	#maybe i dropped singletons in run-mothur?
	
	tmp.table <- count.table[count.table[,i]>0,i,drop=F]
	tmp.table <- tmp.table[rownames(silva97.map.sph)[sample.select.x],] #not sure why this work, rownames have trailing "_"
	map.out <- cbind(silva97.map.sph[sample.select.x,1:2],tmp.table)
	map.out <- cbind(rownames(map.out),map.out)
	colnames(map.out) <- c("name","longitude","latitude","elevation")
	write.table(file=paste0("sample.out.sqrt_edge.",i,".csv"),map.out,sep="\t",quote=F,row.names=F)
}

#export files containing just these samples?


###


map.deg <- cbind(silva97.map.sph,radius.earth+map.alt*100000)
map.xyz <- sph2car(map.deg,deg=T)
map.del <- delaunayn(map.deg)

open3d(windowRect=c(0,0,800,800))
bg3d(sphere=F, back = "filled", color = "lightgrey")
spheres3d(0,0,0, radius = 1.01*radius.earth, color = "lightblue", alpha = 1, back = "fill", front = "fill", lit=F, specular="black", shininess=0, ax.grid = TRUE, sphere.rgl = TRUE)
rgl.sphgrid(radius = 1.05*radius.earth, col.long='red', col.lat='blue', deggap = 15, longtype = "D", add = TRUE, radaxis=FALSE)
map.colors <- rev(terrain.colors(round(max(map.deg[,3])-min(map.deg[,3]))))
map.colors <- map.colors[round(max(map.deg[,3])-map.deg[,3])]
rgl.sphtext(map.deg[,1],map.deg[,2],1.05*map.deg[,3],deg=T,text=".",col=map.colors,cex=1)

mesh.colors <- terrain.colors(round(max(map.deg[,3])-min(map.deg[,3])))
mesh.colors <- mesh.colors[round(max(map.deg[,3]) - apply(matrix(map.deg[map.del,3],ncol=4,byrow=T),1,mean))]
tetramesh(map.del,map.xyz,alpha=1,clear=F,col=mesh.colors,specular="black")



#i <- nndist(silva97.sph[,1:2])
#min(sqrt(i)) #0.17, use grid distance of 0.1
#seq.grid	<- matrix(data=0,nrow=360/0.1,ncol=360/0.1)






