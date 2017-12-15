# Calculate SPEI from climate data

require(SPEI); require(raster); require(rts); require(doParallel); require(foreach); require(iterators); require(ncdf4)

#setwd('~/Documents/TNC/TNC/Data/Statewide/Raster/mortality/california-mortality')
setwd('/Volumes/DGE/Data/Shared/Labs/Field/Private/Marvin/TNC/Statewide/california-mortality')


## General Functions 

	# get latitude raster fxn
	lat_raster = function(x){
		coords = as.data.frame(coordinates(x))
		coords$z = coords$y
		rasterFromXYZ(coords)
	}
	
	# extract dates from raster (LOCA only)
	raster_date_list = function(ras){sapply(names(ras), function(x){
 		as.Date(paste(unlist(strsplit(gsub('X','',x),split="[.]")),collapse='-'))},USE.NAMES = FALSE)}

	
#####################################################	
## LOCA (daily maximum temperature, minimum temperature, and precipitation)
#####################################################
		
# (OPTIONAL) process netcdf data to monthly rasters
		
	# read in data
		had_ppt_85 = stack("Climate/LOCA/hadgem-85-loca-future_ppt.nc")
		had_tmax_85 = stack("Climate/LOCA/hadgem-85-loca-future_tasmax.nc")
		had_tmin_85 = stack("Climate/LOCA/hadgem-85-loca-future_tasmin.nc")
	
	# get monthly sum from daily precip
		had_ppt_85_ts = rts(had_ppt_85,as.Date(raster_date_list(had_ppt_85)))		
		had_ppt_85_monthly = period.apply(had_ppt_85_ts,endpoints(had_ppt_85_ts,'months'),sum)
		
		months = attributes(had_ppt_85_monthly)$time
		monthly_time_series = as.Date(index(months))
		write.csv(monthly_time_series,file="Climate/LOCA/monthly_time_series.csv",row.names=FALSE)
		
		writeRaster(had_ppt_85_monthly@raster, filename="Climate/LOCA/hadgem-85-loca-future_monthly_ppt.tif", 
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))
	
	# get monthly mean from daily min/max		
		had_tmax_85_ts = rts(had_tmax_85,as.Date(raster_date_list(had_tmax_85)))
		had_tmax_85_monthly = period.apply(had_tmax_85_ts,endpoints(had_tmax_85_ts,'months'),mean)
		writeRaster(had_tmax_85_monthly@raster, filename="Climate/LOCA/hadgem-85-loca-future_monthly_tmax.tif", 
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))

		had_tmin_85_ts = rts(had_tmin_85,as.Date(raster_date_list(had_tmin_85)))
		had_tmin_85_monthly = period.apply(had_tmin_85_ts,endpoints(had_tmin_85_ts,'months'),mean)
		writeRaster(had_tmin_85_monthly@raster, filename="Climate/LOCA/hadgem-85-loca-future_monthly_tmin.tif", 
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))
		
		had_tmean_85_monthly = stack(lapply(1:nlayers(had_tmin_85_monthly@raster),function(x)
								calc(stack(had_tmax_85_monthly[[x]], had_tmin_85_monthly[[x]]),mean)))

		writeRaster(had_tmean_85_monthly, filename="Climate/LOCA/hadgem-85-loca-future_monthly_tmean.tif",
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))


# Read in monthly precip and mean temp 	
		had_tmean_85_monthly = stack('Climate/LOCA/hadgem-85-loca-future_monthly_tmean.tif')
		had_ppt_85_monthly = stack('Climate/LOCA/hadgem-85-loca-future_monthly_ppt.tif')
		
		monthly_time_series = as.vector(read.csv("Climate/LOCA/monthly_time_series.csv")[,1])
		
		names(had_tmean_85_monthly) = monthly_time_series
		names(had_ppt_85_monthly) = monthly_time_series
		
# Calculate PET using Thornwaite
	
	# create latitude raster
		raster_lat = lat_raster(had_tmean_85_monthly[[1]])
	
	# parallelize the PET calc over entire timeseries
		cl = makeCluster(detectCores())
		registerDoParallel(cl)
		
		pet_monthly  = stack(
	 					 foreach(i=icount(nlayers(had_tmean_85_monthly)),.packages=c('raster','SPEI'))%dopar%
	 					   calc(stack(had_tmean_85_monthly[[i]], raster_lat),function(y) thornthwaite(y[1],y[2],na.rm=TRUE)))
	
		writeRaster(pet_monthly, filename = "Climate/LOCA/hadgem-85-loca-future_monthly_pet.tif",
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"), overwrite=TRUE)
		
		stopCluster(cl)
	
# Calculate SPEI
		# calc water balance (D)
		spei_d = had_ppt_85_monthly - pet_monthly
		writeRaster(spei_d, filename = "Climate/LOCA/hadgem-85-loca-future_monthly_speiD.tif",
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"), overwrite=TRUE)
		
		# convert to raster timeseries (for viewing)
		spei_d_ts = rts(spei_d,as.Date(monthly_time_series))
		
		# spei over specified time window
		spei_36mo = calc(spei_d_ts@raster, function(x) spei(x,36,na.rm=TRUE)$fitted)
	
		names(spei_36mo) = paste(monthly_time_series)
		
		writeRaster(spei_36mo, filename = "Climate/LOCA/hadgem-85-loca-future_monthly_spei_36mo.tif",
				format="GTiff",datatype='FLT4S',	options=c("COMPRESS=LZW"),overwrite=TRUE)


#####################################################	
# PRISM (monthly total precip, monthly mean temp)
#####################################################

	# PRISM data from earth engine script
		prism_ppt = stack("Climate/PRISM/prism_CA_4km_1990-2017_monthly_ppt_mm.tif")
		prism_tmean_K = stack("Climate/PRISM/prism_CA_4km_1990-2017_monthly_tmean_K.tif")
	
	# convert tmean from deg K to deg C
		prism_tmean = prism_tmean_K - 273 # raster is integer
	
# Calculate PET using Thornwaite

	# get timeseries dates
		yrs = rep(c(1990:2017),each=12)
		mon = rep(1:12,length(c(1990:2017)))
		mon = sprintf("%02d",mon)	
		date_list = paste0(yrs,mon)[-c(323,324,336)] # 3 missing months (11/2016, 12/2016, 12/2017)
	 	monthly_time_series = as.Date(paste0(date_list,"01"),"%Y%m%d")
	 	write.csv(monthly_time_series,file="Climate/PRISM/monthly_time_series.csv",row.names=FALSE)
	
	 
	# convert to raster timeseries (for viewing)
		prism_ppt_ts = rts(prism_ppt, as.Date(paste0(date_list,"01"),"%Y%m%d"))
		prism_tmean_ts = rts(prism_tmean, as.Date(paste0(date_list,"01"),"%Y%m%d"))
	
	# create latitude raster	 
		raster_lat = lat_raster(prism_ppt_ts[[1]])
	
	# parallelize the PET calc over entire timeseries
		cl = makeCluster(detectCores())
		registerDoParallel(cl)
		
		pet_monthly  = stack(
	 					 foreach(i=icount(nlayers(prism_tmean)),.packages=c('raster','SPEI'))%dopar%
	 					   calc(stack(prism_tmean[[i]], raster_lat),function(y) thornthwaite(y[1],y[2],na.rm=TRUE)))
	
		writeRaster(pet_monthly, filename = "Climate/PRISM/prism_4km_1990-2017_monthly_pet.tif",
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))
		
		stopCluster(cl)

# Calculate SPEI for PRISM data
 	# calc water balance (D)
		spei_d = prism_ppt - pet_monthly
		
		writeRaster(spei_d, filename = "Climate/PRISM/prism_4km_1990-2017_monthly_speiD.tif",
				format="GTiff",datatype='FLT4S', options=c("COMPRESS=LZW"))

	# convert to raster timeseries (for viewing)
		spei_d_ts = rts(spei_d,as.Date(paste0(date_list,"01"),"%Y%m%d"))

	# spei over different windows
		#spei_1mo = calc(spei_d, function(x) spei(x,1,na.rm=TRUE)$fitted)
		#spei_12mo = calc(spei_d, function(x) spei(x,12,na.rm=TRUE)$fitted)
		spei_36mo = calc(spei_d, function(x) spei(x,36,na.rm=TRUE)$fitted)
	
		writeRaster(spei_36mo, filename = "Climate/PRISM/prism_4km_1990-2017_monthly_spei_36mo.tif",
				format="GTiff",datatype='FLT4S',	options=c("COMPRESS=LZW"))
	
