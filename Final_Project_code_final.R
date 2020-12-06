#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(sp)
library(maptools)
library(plyr)
library(lubridate)


#Set working directory
dir <- "C:/Users/cam/Documents/GEOG_418/FinalProject/Cameron Anderson/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR("C:/Users/cam/Documents/GEOG_418/FinalProject/Cameron Anderson/Census","Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))


#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 

#Read in dissemination tract shapefile:
census.tracts <- readOGR("C:/Users/cam/Documents/GEOG_418/FinalProject/Cameron Anderson/Census","BC_DA") 
census.tracts <- spTransform(census.tracts, CRS("+init=epsg:26910"))

#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "RdBu", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Income

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=20000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)

####################################################################################
#Descriptive Stats: From lab 1

####For Income####
meanPop <- mean(income.tracts$Income, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation

#Standard Deviation
sdPop <- sd(income.tracts$Income, na.rm = TRUE) #Calculate the SD, ignoring NA values

#Mode
modePop <- as.numeric(names(sort(table(income.tracts$Income), decreasing = TRUE))[1]) #make frequency table of fire size variable and sort it in desending order and extract the first row (Most Frequent)

#Median
medPop <- median(income.tracts$Income, na.rm = TRUE)


#Skewness
skewPop <- skewness(income.tracts$Income, na.rm = TRUE)[1]

#Kurtosis
kurtPop <- kurtosis(income.tracts$Income, na.rm = TRUE)[1]


#CoV
CoVPop <- (sdPop / meanPop) * 100

#Normal distribution test
normPop_PVAL <- shapiro.test(income.tracts$Income)$p.value


#Create a table of descriptive stats

mean = c(meanPop) #Create an object for the means
sd = c(sdPop) #Create an object for the standard deviations
median = c(medPop) #Create an object for the medians
mode <- c(modePop) #Create an object for the modes
skewness <- c(skewPop) #Create an object for the skewness
kurtosis <- c(kurtPop) #Create an object for the kurtosis
CoV <- c(CoVPop) #Create an object for the CoV
normality <- c(normPop_PVAL) #Create an object for the normality PVALUE

##Check table values for sigfigs?

mean
mean <-round(mean,3)

skewness
skewness <-round(skewness, 3)

kurtosis <-round(kurtosis, 3)


sd <-round(sd,3)

mean
sd
median
mode
skewness
kurtosis
normality

####For PM25####

meanPop <- mean(income.tracts$pm2.5, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation

#Standard Deviation
sdPop <- sd(income.tracts$pm2.5, na.rm = TRUE) #Calculate the SD, ignoring NA values

#Mode
modePop <- as.numeric(names(sort(table(income.tracts$pm2.5), decreasing = TRUE))[1]) #make frequency table of fire size variable and sort it in desending order and extract the first row (Most Frequent)

#Median
medPop <- median(income.tracts$pm2.5, na.rm = TRUE)


#Skewness
skewPop <- skewness(income.tracts$pm2.5, na.rm = TRUE)[1]

#Kurtosis
kurtPop <- kurtosis(income.tracts$pm2.5, na.rm = TRUE)[1]


#CoV
CoVPop <- (sdPop / meanPop) * 100

#Normal distribution test
normPop_PVAL <- shapiro.test(income.tracts$pm2.5)$p.value


#Create a table of descriptive stats

mean = c(meanPop) #Create an object for the means
sd = c(sdPop) #Create an object for the standard deviations
median = c(medPop) #Create an object for the medians
mode <- c(modePop) #Create an object for the modes
skewness <- c(skewPop) #Create an object for the skewness
kurtosis <- c(kurtPop) #Create an object for the kurtosis
CoV <- c(CoVPop) #Create an object for the CoV
normality <- c(normPop_PVAL) #Create an object for the normality PVALUE

##Check table values for sigfigs?

mean
mean <-round(mean,3)

skewness
skewness <-round(skewness, 3)

kurtosis <-round(kurtosis, 3)


sd <-round(sd,3)

mean
sd
median
mode
skewness
kurtosis
normality

###################################################################################
#Global/Local Moran's I: Obj 1, spatial segregation of income.(lab 3)

inc.i <- poly2nb(income.tracts)               #CALCULATES NEIGBOUR WEIGHTS MATRIXS. POLYGONS THAT ARE NEIGHBOURS ARE VALUED 1
crd.net <- nb2lines(inc.i,coords=coordinates(income.tracts)) #CONVERTS NEIGHBOUR MATRIX INTO LINES WE CAN PLOT


inc.lw <- nb2listw(inc.i, zero.policy = TRUE, style = "W")#WEIGHT MATRIX
print.listw(inc.lw, zero.policy = TRUE)

mi <- moran.test(income.tracts$Income, inc.lw, zero.policy = TRUE)#LOOKS AT MEDIAN INC. GIVE IT WIEGHT MATRIX
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(inc.lw)

mI <- mi$estimate[[1]] #MORANS I
eI <- mi$estimate[[2]] #EXPECTED
var <- mi$estimate[[3]] #VARIANCE

z <- ((mI-(eI))/sqrt(var))#IS THE MORANS I SIGNIFICANT
mI
eI
z
########################  
lisa.test <- localmoran(income.tracts$Income, inc.lw)#GIVE IN MEDIAN INC AND WEIGHT MATRIX

income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

map_LISA <- tm_shape(income.tracts) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fisher", 
              palette = "-RdBu", n = 6, midpoint=NA) +
  tm_legend(legend.outside=TRUE)


map_LISA

moran.plot(income.tracts$Income, inc.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Population Density", 
           ylab="Spatially Lagged Population Density", quiet=NULL)



##################################################################################
#Spatial Interpolation (lab 4)


proj4string(grd) <- proj4string(pm2.5)#use same projection
P.idw <- gstat::idw(PM25 ~ 1, pm2.5, newdata=grd, idp=2)
r      <- raster(P.idw)                                             
r.m     <- mask(r, census.tracts)



tm_shape(r.m) + 
  tm_raster(n=10,palette = "-RdBu",
            title="Predicted PM 2.5 \n(in ppm)") + 
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



# Leave-one-out validation routine ###check how it woks
IDW.out <- vector(length = length(pm2.5))
for (i in 1:length(pm2.5)) {                                                       
  IDW.out[i] <- gstat:: idw(PM25 ~ 1, pm2.5[-i,], pm2.5[i,], idp=2)$var1.pred
}#leaves out a point, fit a surface, compare it to original. )

OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ pm2.5$PM25, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ pm2.5$PM25), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt( sum((IDW.out - pm2.5$PM25)^2) / length(pm2.5))#RSME 

############################

# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(PM25~1, pm2.5, newdata=grd, idp=2)
n   <- length(pm2.5)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(PM25~1, pm2.5[-i,], newdata=grd, idp=2)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 


r <- raster(img.sig, layer="v")
r.m <- mask(r, census.tracts)

#Plot the map
tm_shape(r.m) +
  tm_raster(style = "fixed",
            auto.palette.mapping = FALSE,
            palette = "YlOrRd", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, Inf),
            title="95% confidence interval \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

#################################################################################

#Point Pattern (lab 2)

###NEAREST NEIGHBOUR
kma <- pm2.5


kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]
#check for and remove duplicated points
#check for duplicated points
#finds zero distance among points
zd <- zerodist(kma)
zd
#remove duplicates
kma <- remove.duplicates(kma)

#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(kma)) #stick the extent into a matrix 
kma.ext

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
View(kma.ppp)

census.tracts <- spTransform(census.tracts, CRS("+init=epsg:3005"))

##Nearest Neighbour Distance
###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(kma.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
View(pm2.5)
nnd = (sum(nearestNeighbour$Distance))/306

#mean nearest neighbour for random spatial distribution

studyArea <- area(census.tracts)
pointDensity <- 306/studyArea

r.nnd = 1/(2*sqrt(pointDensity))

d.nnd = 1.07453/(sqrt(pointDensity))

R = nnd/r.nnd

SE.NND <- 0.26136/(sqrt(306*pointDensity))

z = (nnd-r.nnd)/SE.NND
R
z
r.nnd
d.nnd
nnd

################################################################################
#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
# Convert your interpolation into a raster and map it:
r <- raster(r.m)
sufaceMap <- tm_shape(r.m) + 
  tm_raster(n=5,palette = "RdBu",
            title="PM 2.5 \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2)
sufaceMap
#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average pm2.5 for each polygon
income.tracts$pm2.5 <- round(extract(r.m, income.tracts, fun=mean)[,1], 5)

################################################################################
######Linear Regression##########
#Let's say your dataset with both PM2.5 and Income 
#are stored in a dataset called income.tracts.
#Plot income and PM2.5 from the income.tracts dataset you created
plot(income.tracts$Income~income.tracts$pm2.5)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts.no0 <-  income.tracts[which(income.tracts$pm2.5 > 0), ]

#Now plot the data again
plot(income.tracts.no0$Income~income.tracts.no0$pm2.5,
     xlim=c(0,0.5))

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$pm2.5,
               xlim=c(0,0.5))
#Add the regression model to the plot you created
#Graph worked best first time so use that one
plot(income.tracts$Income~income.tracts$pm2.5,
     xlim=c(0,0.5))
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts.no0)

#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Residuals",
              style = "jenks",
              palette = "-RdBu", n = 6) +
  tm_legend(legend.outside=TRUE)

map_resid

#Global Moran's I

resids.i <- poly2nb(income.tracts.no0)                              #CALCULATES NEIGBOUR WEIGHTS MATRIXS. POLYGONS THAT ARE NEIGHBOURS ARE VALUED 1
crd.net <- nb2lines(inc.i,coords=coordinates(income.tracts)) #CONVERTS NEIGHBOUR MATRIX INTO LINES WE CAN PLOT


resids.lw <- nb2listw(resids.i, zero.policy = TRUE, style = "W")#WEIGHT MATRIX
print.listw(inc.lw, zero.policy = TRUE)

mi <- moran.test(income.tracts.no0$residuals, resids.lw, zero.policy = TRUE)#LOOKS AT MEDIAN INC. GIVE IT WIEGHT MATRIX
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(resids.lw)

mI.r <- mi$estimate[[1]]#MORANS I
eI.r <- mi$estimate[[2]]#EXPECTED
var.r <- mi$estimate[[3]]#VARIANCE

z.r <- ((mI.r-(eI.r))/sqrt(var.r))#IS THE MORANS I SIGNIFICANT

mI.r
eI.r
z.r
 

lisa.test.r <- localmoran(income.tracts.no0$residuals, resids.lw)#GIVE IN MEDIAN INC AND WEIGHT MATRIX

income.tracts.no0$Ii <- lisa.test.r[,1]
income.tracts.no0$E.Ii<- lisa.test.r[,2]
income.tracts.no0$Var.Ii<- lisa.test.r[,3]
income.tracts.no0$Z.Ii<- lisa.test.r[,4]
income.tracts.no0$P<- lisa.test.r[,5]

resids.lw <- nb2listw(resids.i, zero.policy = TRUE, style = "W")#WEIGHT MATRIX


map_LISA <- tm_shape(income.tracts.no0) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fisher", 
              palette = "-RdBu", n = 5,
              midpoint= NA) +
  tm_legend(legend.outside=TRUE)


map_LISA

moran.plot(income.tracts$Income, inc.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Residuals", 
           ylab="Spatially Lagged Residuals", quiet=NULL)





#############################################################################
####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)
#Observe the result:
head(income.tracts.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "-RdBu", n = 6, midpoint= NA) +
  tm_legend(legend.outside=TRUE)
map_r2

#Time for more magic. Let's map the coefficients
income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5
#Create choropleth map of the coefficients
map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "-RdBu", n = 6, midpoint= NA) +
  tm_legend(legend.outside=TRUE)
map_coef






