objects()
rm(list=objects())
#
# Build a working directory for this weeks lab and change working dir
# Note you might have to specify the path explicitly  as some 
# computers in the lab were not working correctly, to do this go to
# >Session>Set Working Directory
dir.create("~/Week07")
setwd("~/Week07")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr)

# Get our gold standard flow data from USGS 0205551460 `
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                           end_date = "2019-01-01")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# We can determine our UTM zone and build our proj4 projection string
trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)

myflowgage$area   # area in km2
[1] 12.9792
# If the watershed was square, which it is not, the size would be the 
# square root of the area. Also, the gage/pour point is not in the center
# so we will search around the gage.
# Build sp point for USGS gage location, in both _ll and _utm
latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
             myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)
# We are going to over estimate our area
sqrt(myflowgage$area)   # guestimating square watershed
# For our search we are going to multiply the area by 6 and
# to get the distance
sqrt(myflowgage$area*8)
searchlength=sqrt(myflowgage$area*8)*1000 
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")




#zoom(mydem)
#plot(bboxpts,add=T)
#plot(pourpoint,add=T,col="red")

old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                        sep = ":"))
system("mpirun aread8")
# Write our raster to a geotiff file that can be used with
# OS level hydrological models 
writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd;  #Linux/Mac
# dir; cd ; 

z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)

plot(z-fel)
# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
p=raster("mydemp.tif")
plot(p)
sd8=raster("mydemsd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 8 aread8 -p mydemp.tif -ad8 mydemad8.tif")

ad8=raster("mydemad8.tif")
plot(log(ad8))
#zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
plot(ang)
slp=raster("mydemslp.tif")
plot(slp)


# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
plot(log(sca))
#zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 2000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint,add=T)
zoom(src)
# a quick R function to write a shapefile  #how much water is going
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

zoom(src)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 4000")
src1=raster("mydemsrc1.tif")
plot(src1)
zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"),add=T)
#zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"),add=T)

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)
#zoom(ssa)
#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
zoom(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))



##Copyying from lab notes
# Now that we can get Flow Accumulation clipped to the DEM
# and as in previous labs we can get slope from the clipped DEM
# we can use our lecture notes to calculate the Topographic Index(TI)
# 
# Trimming, Cropping, and Masking to make life prettier and easier

mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
zoom(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
#For difference between DEM and filled DEM
plot(z-fel)
mybasindif=crop(fel,mybasinmask)
mybasindif=mask(mybasindif,mybasinmask)
plot(mybasindif-mybasindem)
# A series of plots to show all of the components
#
par(mfrow = c(2, 2))
plot(TI)
plot(TIC)
plot(mybasindem)
plot(mybasindif-mybasindem)
#plot(mybasinsca)
#plot(mybasinslp)
dev.off()
plot(TIC,col=rainbow(5))
# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T)
plot(mydemw_poly,add=T,border="red")
mydemw_poly
#writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab
# Open the WebSoilSurvey site to: 
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/32bhl1fdbpn253bd0sln1pkz/wss_aoi_2022-03-17_11-40-17.zip"
download.file(url,"wss_aoi.zip")
unzip("wss_aoi.zip")
mysoil=readOGR("wss_aoi_2022-03-17_11-40-17/spatial/soilmu_a_aoi.shp") 

list.files()
list.files(pattern = "wss")
list.files("wss_aoi_2022-03-17_11-40-17/spatial/",pattern = "shp")
list.files(paste0(list.files(pattern="wss"),"/spatial/"))

# This needs to be completed based on your download
mysoil=readOGR("wss_aoi_2022-03-17_11-40-17/spatial/soilmu_a_aoi.shp")    
# Explore the mysoil dataset which is returned
head(mysoil@data)
class(mysoil)
mysoil
dev.off()
plot(mybasindem)
#proj4string(mysoil) = proj4_ll
#mysoil1=spTransform(mysoil,crs_utm)
plot(mysoil1)
# https://www.usgs.gov/centers/eros/science/national-land-cover-database

#other way
library(rgdal)
library(raster)
library(rgeos)
#system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
#mydemw_poly_gdal=readOGR("mydemw_poly_gdal.dbf")
#x <- mapunit_geom_by_ll_bbox(bbox(spTransform(mydemw_poly,crs_ll)))
#x=spTransform(x,crs_utm) # convert LL to UTM
#out <- gIntersection(x, mydemw_poly, byid=TRUE)
#plot(out,col=rainbow(100))

par(mfrow = c(2, 2))
plot(TI, main="TI")
plot(TIC, main="TIC")
plot(mybasindem,main="DEM")
plot(mysoil,main="Soils",col=rainbow(100))

#install.packages("FedData")
library(FedData)
#Adding landcover
vepPolygon <- polygon_from_extent(raster::extent(588297.3, 594362, 4126215, 4133068), 
                                  proj4string= '+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs')                                                                                       
NLCD <- get_nlcd(template=vepPolygon, label='VEPIIN')
