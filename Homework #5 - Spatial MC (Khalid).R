# III. INTRADISTRIBUTIONAL MOBILITY WITH MARKOV CHAINS (MC) ####

install.packages('estdaR') # install package from CRAN
install.packages('ggmap')
install.packages('ggplot2')
install.packages('maptools')
install.packages('spdep')
install.packages('devtools') # collection of package development tools
devtools::install_github("amvallone/estdaR", force=TRUE) # force install. from Github

library(rgdal) 
library(ggmap)
library(ggplot2)
library(maptools)
library(spdep)
library(devtools)
library(estdaR)
library(sf)
library(sp)

# III.2. Read and manage database ####

nigeria <- read.csv("D:/UNIVERSITY 2.04/UBFC Dijon 2.04/UBFC study/semester 2/4. Advanced Topics in Spatial Statistics/nigeria.population.csv", header=TRUE)

b00 <- as.data.frame(nigeria) # Coerce a vector map into double entry table
b0 <- as.data.frame(b00[,4:14]) # Select years 1990-2018
b <- apply(b0,2,as.numeric) # Coerce cols. into numeric (integer/real)

# statenames<-c("1QT","2QT","3QT","4QT","5QT") # forms a row vector
statenames<-c("very small","small","average","high","highest") # forms a row vector

set.seed(2010)

# III.3. Plot a map ####

coordinates(nigeria) <- ~latitude+longitude
cor <- coordinates(nigeria)

lon <- nigeria$longitude
lat <- nigeria$latitude

bbox <- make_bbox(lon,lat,f=0.05)
map.nigeria <- get_map(bbox,maptype = "terrain")

m1 <- ggmap(map.nigeria) 
m1

# III.4. Spatial weight matrix, W ####

# Nearest neighbor W matrix
k1 <- knn2nb(knearneigh(cor, k=1)) # spdep package
k1d <- nbdists(k1, cor)
max_k1 <- max(unlist(k1d))
dnn_nb <- dnearneigh(cor, 0, max_k1)
summary(dnn_nb)
plot(nigeria)
plot(dnn_nb, cor, add=TRUE)

# Inverse distance nearest neighbor W matrix
deuc_dnn<-nbdists(dnn_nb,cor)
dinv_dnn <- lapply(deuc_dnn, function(x) 1/x)

# 3-nearest neighbors W matrix
knn3_nb <- knn2nb(knearneigh(cor, k=3),
                  row.names=row.names(nigeria))
summary(knn3_nb)
plot(nigeria)
plot(knn3_nb, cor, add=TRUE)

# Standardized spatial weight matrices, W*
dmins <-nb2listw(dnn_nb,style="W")
dinvs <- nb2listw(dnn_nb, glist=dinv_dnn)
knn3s <-nb2listw(knn3_nb,style="W")

# III.5. Spatial Markov Transition matrix ####

spatial.markov <- sp.mkv(b,dinvs) # spatial Markov chain
p<-spatial.markov[[1]] # probabilities
t<-spatial.markov[[2]] # transitions
p
t
St <- do.call(rbind, lapply(1:5,function(x)st.st(spatial.markov[[1]][,,x])))
St

# Table with heatmap
# Prepare a full matrix binding each matrix with its ergodic
m<-rbind(p[,,1],St[1,],p[,,2],St[2,],p[,,3],St[3,],p[,,4],St[4,],p[,,5],St[5,])
# create a set of vector to improve the heat map plot
l.l <- c(rep("l(1) ",6),rep("l(2) ",6),rep("l(3) ",6),rep("l(4) ",6),rep("l(5) ",6))
rownames(m) <- paste(l.l, rep(c(paste("Q ",1:5,sep=""), "\u03C0"),5),sep=" ")
colnames(m) <- paste("Q",1:5,sep=" ") #create a label vector for the heat map
col<-c("white","red") #set the colors
m<-round(m,3) 
data <- reshape2::melt(m) # reshape the matrix to a long version data frame 
names(data) <- c("y","x","Probabilities")
s <- rep(.4,150) #adjust the text size
s[seq(6,150,6)]<-.7

h.cities <- ggplot(data,aes(x,y,Probabilities,label=sprintf("%0.3f", round(Probabilities, digits = 2)))) +
  geom_tile(aes(fill=Probabilities),color="black",size=s)+theme_bw() +
  scale_fill_gradient(low=col[1], high=col[2]) +
  labs(x=NULL, y=NULL,title=NULL) +
  theme(axis.title=element_text(size=14,face="bold"),
        panel.border = element_rect(colour="white"),
        axis.text = element_text(size=12),
        axis.line = element_line(),
        legend.text =element_text(size=12),
        legend.title=element_text(size=.8*14,face="bold"),
        plot.title = element_text(size=2*14,hjust=0.5))+
  geom_text(size=4)+scale_x_discrete(position = "top")+
  scale_y_discrete(limits = rev(levels(data$y)))
h.cities

# Ergodic distributions
p1 <- sp.mkv(b,classes=5,dinvs)[[1]] #idem a mkv pero condicionada espacialmente.
Sergo <- do.call(rbind, lapply(1:5,function(x)st.st(p1[,,x]))) #computes the ergodic distribution for each conditional matrix 
colnames(Sergo)<-statenames # Set the column names of a matrix-like object.
splags<-c("l(1)","l(2)","l(3)","l(4)","l(5)")
rownames(Sergo)<-splags
Sergo


# Test for homogeneity of the spatial Markov transition matrix
smct<-sp.homo.test(b,dinvs,fixed=FALSE)
smct
nullm<-smct[[3]]
colnames(nullm)<-statenames
rownames(nullm)<-statenames
nullm

# III.6. LISA Markov transition matrix ####

L.M <- lisamkv(b,dinvs,geoda=FALSE)
mLM<-L.M[["move"]] # move type from 1 to 16
tLM<-L.M[["lisamatrix"]] # transition matrix
pLM<-L.M[["p.lisamatrix"]] # probability matrix
lisacat<-c("HH","LH","LL","HL")
colnames(tLM)<-lisacat
rownames(tLM)<-lisacat
tLM

colnames(pLM)<-lisacat
rownames(pLM)<-lisacat
pLM

gro<-c("2006","2007","2008","2009","2010","2011","2012","2013","2014","2015")
colnames(mLM)<-gro
mLM

# Ergodic distribution
lisaer<-st.st(L.M$p.lisamatrix)
lisaer

# Test of co-movement dependence - provisional
# cmov<-join.d(b,dinvs)
# mcmov<-cmov[["Expected"]] # LISA markov transition matrix under H0
# colnames(mcmov)<-lisacat
# rownames(mcmov)<-lisacat
# mcmov

# III.7. Directional LISA ####

# Spatial regimes: Center, North, South
Regime <- rep("East", dim(nigeria)[1]) # Replic. "South" in row dim. of 'nigeria'. 
Regime[which(nigeria$FIELD2=="I"| nigeria$FIELD2=="XV"| nigeria$FIELD2=="II" | nigeria$FIELD2=="III" | nigeria$FIELD2 =="IV")] <- "North"

# #log. transformation and mean normalization
den<-function(x){
  y<-log(x)
  z<-mean(y)
  out<-y/z
  out
}

D1 <- apply(b,2,den)
D <- D1

# Moran scatterplot and rose diagram
full.p8 <- d.LISA(D[,1],D[,8],dinvs,Regime=Regime,nsim=999)
full.p4 <- d.LISA(D[,1],D[,8],dinvs,Regime=Regime,nsim=999,k=4)
full.p16 <- d.LISA(D[,1],D[,8],dinvs,Regime=Regime,nsim=999,k=16)