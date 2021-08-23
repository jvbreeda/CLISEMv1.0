###############################################################################################################
# main routine
###############################################################################################################

#load netcdf libraries and source PCA emulator code
library(ncdf4)
source('pca_emul.R')

#to normalize the input data (saved during the calibration process)
load("sx1.Rda")
load("sx2.Rda")
load("sx3.Rda")
load("sx4.Rda")
load("sx5.Rda")

#----------------------------------------------------------------------------------------------------
#	READ THE INPUT DATA: orbital parameters, CO2 and ice sheet parameter
#----------------------------------------------------------------------------------------------------
x_full=read.table('emul_input.txt') #read the input data (orbital parameters and CO2 values)
x_vol=scan('EMULICE')  		    #read the input data: ice sheet paramater (ice volume here)
index=x_vol[[1]]		    #index of the year at which we would extract the variables from the matrix x_full
x_ice=x_vol[[2]]    		    #ice volume as the fifth variable

x_year=x_full[index,]		    #extract the proper year

var1=(x_year$V1)/sx1		    #read the value of esinw at the current year and normalize
var2=(x_year$V2)/sx2		    #read the value of ecosw at the current year and normalize
var3=(x_year$V3)/sx3		    #read the value of obliquity at the current year and normalize
var4=(x_year$V4)/sx4		    #read the value of CO2 at the current year and normalize
var5=x_ice/sx5			    #read the value of the ice sheet parameter (ice volume here) at the current year and normalize

# make vector of all input variables at the current year
x=cbind(var1,var2,var3,var4,var5)



#----------------------------------------------------------------------------------------------------
#	TEMPERATURE PREDICTION
#----------------------------------------------------------------------------------------------------
#load the PC scores of the calibration for the temperatures

prefix="E_list_temp"
suffix=".Rda"


load(paste(prefix,suffix,sep=""))

# estimate the current climate using the PCA routine pe_p 
estimate_temperature <- function (x,  E_list ) 
{
  lapply (  E_list ,  function(E)   pe_p (x, E) )
}

temperature <- estimate_temperature(x, E_list_temp)

# the mean temperature field for each month is:
temp_jan = temperature[[1]]$mean-273.15 
temp_dec = temperature[[12]]$mean-273.15 

# create an array on the ice sheet model grid of the form (LAT, LON) : 
yval=201
xval=201
nmonths=12

OUTPUT_TEMP <- array ( temp_jan, dim  = c(yval, xval, nmonths))

for (imonth in seq(1,12))
{
  OUTPUT_TEMP [ , ,imonth] = temperature[[imonth]]$mean-273.15
}


#----------------------------------------------------------------------------------------------------
#	PRECIPITATION PREDICTION
#----------------------------------------------------------------------------------------------------

prefix="E_list_precip"
suffix=".Rda"
load(paste(prefix,suffix,sep=""))

# estimate the current climate using the PCA routine pe_p 
estimate_precipitation <- function (x,  E_list ) 
{
  lapply (  E_list ,  function(E)   pe_p (x, E) )
}

precipitation <- estimate_precipitation(x, E_list_precip)

# the mean precipitation field is: 86400 s per day, 30 days per month (from kg m-2 s-1 = mm s-1 to m of water equivalent per month)
prec_jan= precipitation[[1]]$mean*86400*30/1000
prec_dec= precipitation[[12]]$mean*86400*30/1000

# create an array on the ice sheet model grid of the form (LAT, LON) : 
yval=201
xval=201

OUTPUT_PREC <- array ( temp_jan, dim  = c(yval, xval, nmonths))

for (imonth in seq(1,12))
{
  OUTPUT_PREC [ , ,imonth] = precipitation[[imonth]]$mean*86400*30/1000
}


#----------------------------------------------------------------------------------------------------
# write temperature to NetCDF file
#----------------------------------------------------------------------------------------------------
name="MONTHLY_TEMP"
nc=".nc"

xval = ncdim_def("lon","degrees_north",seq(1:201))
yval = ncdim_def("lat","degrees_north",seq(1:201))
dims_exp = ncdim_def("months","number",seq(1:12))
temp_all = ncvar_def("monthly T","degrees K",list(xval,yval,dims_exp))
ncnew    = nc_create(paste(name,nc,sep=""),temp_all)
ncvar_put(ncnew,temp_all,OUTPUT_TEMP)
nc_close(ncnew)

name="MONTHLY_PRECIP"
nc=".nc"

xval = ncdim_def("lon","degrees_north",seq(1:201))
yval = ncdim_def("lat","degrees_north",seq(1:201))
dims_exp = ncdim_def("months","number",seq(1:12))
precip_all = ncvar_def("monthly precipitation","m rainfall",list(xval,yval,dims_exp))
ncnew    = nc_create(paste(name,nc,sep=""),precip_all)
ncvar_put(ncnew,precip_all,OUTPUT_PREC)
nc_close(ncnew)




