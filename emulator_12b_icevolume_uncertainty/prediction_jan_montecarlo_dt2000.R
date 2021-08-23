##################################################################################################################################
# STOCHASTIC MODEL USED FOR THE MONTE CARLO SIMULATIONS
# 
##################################################################################################################################

library(ncdf4)
source('pca_emul.R')
#require(BBmisc)

##################################################################################################################################
#			PREDICTION - DETERMINISTIC MODEL
##################################################################################################################################
load("sx1.Rda")
load("sx2.Rda")
load("sx3.Rda")
load("sx4.Rda")
load("sx5.Rda")

x_full=read.table('emul_input.txt')		 #read the orbital parameters (esinw, ecosw, obliquity) and CO2 concentration
x_vol=scan('EMULICE')			         #read the ice volume and the year at the present time step
x_vol_full=read.table('emulice_full_dt2000.txt')  #read the ice volume and the year for the entire length of the experiment

index=x_vol[[1]]    			     	 #index at the present time step, the indices start at 40 Myr ago and are given every 250 years
						 # since the experiments start at 38 Myr ago, the first index is 8001 (2,000,000/250=8000)
x_ice=x_vol[[2]]    			      	 #ice volume as the second variable from the file EMULICE
x_ice_1=x_vol_full[[(index-8001)/8,2]]      	 #ice volume time step -1, index-8001/2 for dt = 500 years and index-8001/8 for dt = 2000 years


index_1=index-8				      	 # previous index is index-2 for dt = 500 years and index-8 for dt = 2000 years

x_year=x_full[index,]			      	 #extract the proper year
x_year_1=x_full[index_1,]			 #extract the proper year -1


var1=(x_year$V1)/sx1                		 #read the value of esinw at the current time step and normalize
var2=(x_year$V2)/sx2                		 #read the value of ecosw at the current time step and normalize
var3=(x_year$V3)/sx3                		 #read the value of obliquity at the current time step and normalize
var4=(x_year$V4)/sx4                		 #read the value of CO2 at the current time step and normalize
var5=x_ice/sx5                      		 #read the value of the ice sheet parameter (ice volume here) at the current time step and normalize

var1_1=x_year_1$V1/sx1		    		 #read the value of esinw at the previous time step and normalize
var2_1=x_year_1$V2/sx2		    		 #read the value of ecosw at the previous time step and normalize
var3_1=x_year_1$V3/sx3		    		 #read the value of obliquity at the previous time step and normalize
var4_1=x_year_1$V4/sx4		    	         #read the value of CO2 at the previous time step and normalize
var5_1=x_ice_1/sx5		    	         #read the value of the ice sheet parameter (ice volume here) at the previous time step and normalize


x=cbind(var1,var2,var3,var4,var5)		 #input variables at time step i
x_1=cbind(var1_1,var2_1,var3_1,var4_1,var5_1)	 #input variables at time step i-1

x_stoch=matrix(c(var1_1, var2_1, var3_1, var4_1, var5_1, var1, var2, var3, var4, var5), ncol=5, nrow=2,byrow=TRUE)	# input variables at time step i and time step i-1

#load the PCA emulator source scripts
source('./pca_emul.R')

# ----------------------------------------------------------------------------------------------------
#			CALL THE STOCHASTIC MODEL AND CALCULATE THE JANUARY TEMPERATURE
# ----------------------------------------------------------------------------------------------------

# load the output from the previous time step to call the function "stochastic_pca_onestep"
# calibrated with only January temperayures; contains E_List_jan2.Rda
load('E_stochastic.Rda')


#only execute at the first time step of the simulations
if((index-8001)/2 < 2){
stoch1 <- stochastic_pca_onestep_null_oldinput(E_list_jan2,x)

# execute at all other timesteps of the simulations
# E_list_jan2 contains the mean, variance and calibrated PC scores from the calibration process
# el1 contains the simulated PC scores from time step -1
} else {
load('el1.Rda')
#stoch1 <- stochastic_pca_onestep(E_list_jan2,x,el1,x1)
stoch1 <- stochastic_pca_onestep(E_list_jan2,x_stoch,el1)
}

# save the output from the current time step, to call during the next time step
save(stoch1,file="el1.Rda")

# the mean field for january is
data= stoch1$simulated.field-273.15 

# create an array the size of the ice sheet model grid
# the ice sheet model at 40 km resolution has 201 points in x and y 
yval=201
xval=201
OUTPUT <- array ( data, dim  = c(yval, xval))

# ----------------------------------------------------------------------------------------------------
#			write new January temperature to a NetCDF file
# ----------------------------------------------------------------------------------------------------

name="temp"
nc=".nc"
month = "jan"

xval = ncdim_def("lon","degrees_north",seq(1:201))
yval = ncdim_def("lat","degrees_north",seq(1:201))
temp_all = ncvar_def("monthly T","degrees K",list(xval,yval))
ncnew    = nc_create(paste(name,month,nc,sep=""),temp_all)
ncvar_put(ncnew,temp_all,OUTPUT)
nc_close(ncnew)

