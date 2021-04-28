# -------------------------------------------------------------------------------------------------------------------------
# read the temperature and precipitation data from the HadSM3 runs
# -------------------------------------------------------------------------------------------------------------------------
#load libraries
library(ncdf4)
library(fields)
library(maps)

#define the experiments names (100 experiments: xaemaa -> xaemdv)
Exp1="xaema"
Exp2="xaemb"
Exp3="xaemc"
Exp4="xaemd"

Exp_list_A = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
Exp_list_B = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
Exp_list_C = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
Exp_list_D = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v")

month="aug"
prefix="mean_pd_"
suffix=".nc"

# -------------------------------------------------------------------------------------------------------------------------
# READ IN INPUT DATA = orbital parameters,CO2 and the ice sheet parameter (ice volume here)
# -------------------------------------------------------------------------------------------------------------------------
cont_paramdat = read.table('emul_input_5variables.txt', header=TRUE)
cont_param_dim=dim(cont_paramdat)
obl = cont_paramdat$obliquity
esinw = cont_paramdat$esinw
ecosw = cont_paramdat$ecosw
co2 = cont_paramdat$co2_ppm
ice = cont_paramdat$ice



# --------------------------------------------------------------------------------------------------------------------------
# READ IN OUTPUT DATA = temperature and precipitation
# --------------------------------------------------------------------------------------------------------------------------
temp=list((0))
precip=list((0))
# read the latitude and longitude on the ice sheet model grid
       filepathGRID=file.path('..','output_emulator_100a',"AISM40_eot.nc")
       ncgrid = nc_open(filepathGRID,write=TRUE)
       lats_dat = ncvar_get(ncgrid, "lat")
       lons_dat = ncvar_get(ncgrid, "lon")

# read the air temperature and precipitation from the climate model runs 
for (i in Exp_list_A) {
       	filepathA=file.path('..','output_emulator_100a',paste(Exp1,i,sep=""),paste(prefix,month,suffix,sep=""))
	nc = nc_open(filepathA,write=TRUE)
        temp = ncvar_get(nc, "temp")
	nam = paste("temp_xaema", i, "_kdat", sep = "")
	assign(nam, temp)
        precip = ncvar_get(nc, "precip")
	nam_precip = paste("precip_xaema", i, "_kdat", sep = "")
	assign(nam_precip, precip)	
	nc_close(nc)
	rm(nc, temp, precip)
}

for (i in Exp_list_B) {
        filepathB=file.path('..','output_emulator_100a',paste(Exp2,i,sep=""),paste(prefix,month,suffix,sep=""))
    	nc = nc_open(filepathB,write=TRUE)
      	temp = ncvar_get(nc, "temp")
      	nam = paste("temp_xaemb", i, "_kdat", sep = "")
      	assign(nam, temp)
        precip = ncvar_get(nc, "precip")
	nam_precip = paste("precip_xaemb", i, "_kdat", sep = "")
	assign(nam_precip, precip)	
      	nc_close(nc)
      	rm(nc, temp, precip)
}

for (i in Exp_list_C) {
        filepathC=file.path('..','output_emulator_100a',paste(Exp3,i,sep=""),paste(prefix,month,suffix,sep=""))
    	nc = nc_open(filepathC,write=TRUE)
      	temp = ncvar_get(nc, "temp")
      	nam = paste("temp_xaemc", i, "_kdat", sep = "")
      	assign(nam, temp)
        precip = ncvar_get(nc, "precip")
	nam_precip = paste("precip_xaemc", i, "_kdat", sep = "")
	assign(nam_precip, precip)	
      	nc_close(nc)
      	rm(nc, temp, precip)
}

for (i in Exp_list_D) {
       filepathD=file.path('..','output_emulator_100a',paste(Exp4,i,sep=""),paste(prefix,month,suffix,sep=""))
       nc = nc_open(filepathD,write=TRUE)
       temp = ncvar_get(nc, "temp")
       nam = paste("temp_xaemd", i, "_kdat", sep = "")
       assign(nam, temp)
       precip = ncvar_get(nc, "precip")
       nam_precip = paste("precip_xaemd", i, "_kdat", sep = "")
       assign(nam_precip, precip)	
       nc_close(nc)
       rm(nc, temp, precip)
}

# --------------------------------------------------------------
# MAKE MISSING VALUES
# --------------------------------------------------------------

for (i in Exp_list_A) {
  nam = paste("temp_xaema", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)

  nam_precip = paste("precip_xaema", i, "_kdat", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_precip = replace(data_precip, data_precip == 2.00000004008175e+20, NA)
  assign(nam_precip, data_precip)
  rm(data_precip)

}

for (i in Exp_list_B) {
  nam = paste("temp_xaemb", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)

  nam_precip = paste("precip_xaemb", i, "_kdat", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_precip = replace(data_precip, data_precip == 2.00000004008175e+20, NA)
  assign(nam_precip, data_precip)
  rm(data_precip)
}

for (i in Exp_list_C) {
  nam = paste("temp_xaemc", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)

  nam_precip = paste("precip_xaemc", i, "_kdat", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_precip = replace(data_precip, data_precip == 2.00000004008175e+20, NA)
  assign(nam_precip, data_precip)
  rm(data_precip)

}

for (i in Exp_list_D) {
  nam = paste("temp_xaemd", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)

  nam_precip = paste("precip_xaemd", i, "_kdat", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_precip = replace(data_precip, data_precip == 2.00000004008175e+20, NA)
  assign(nam_precip, data_precip)
  rm(data_precip)

}


# ----------------------------------------------------------
# VECTORIZE TEMPERATURE MATRICES (BY COLUMN (LATITUDE))
# ----------------------------------------------------------
temp_tdum_k_all = matrix(0, nrow = length(temp_xaemaa_kdat), ncol = 100)
precip_tdum_k_all = matrix(0, nrow = length(temp_xaemaa_kdat), ncol = 100)


j = 1

for (i in Exp_list_A) {
  nam = paste("temp_xaema", i, "_kdat", sep = "")
  nam_new = paste("temp_xaema", i, "_k_vec", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data_new = as.vector(data)
  assign(nam_new, data_new)
  temp_tdum_k_all[,j] = data_new
  rm(data)

  nam_precip = paste("precip_xaema", i, "_kdat", sep = "")
  nam_new_precip = paste("precip_xaema", i, "_k_vec", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_new_precip = as.vector(data_precip)
  assign(nam_new_precip, data_new_precip)
  precip_tdum_k_all[,j] = data_new_precip
  rm(data_precip)

  j = j + 1
}

for (i in Exp_list_B) {
  nam = paste("temp_xaemb", i, "_kdat", sep = "")
  nam_new = paste("temp_xaemb", i, "_k_vec", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data_new = as.vector(data)
  assign(nam_new, data_new)
  temp_tdum_k_all[,j] = data_new
  rm(data)

  nam_precip = paste("precip_xaemb", i, "_kdat", sep = "")
  nam_new_precip = paste("precip_xaemb", i, "_k_vec", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_new_precip = as.vector(data_precip)
  assign(nam_new_precip, data_new_precip)
  precip_tdum_k_all[,j] = data_new_precip
  rm(data_precip)

  j = j + 1
}


for (i in Exp_list_C) {
  nam = paste("temp_xaemc", i, "_kdat", sep = "")
  nam_new = paste("temp_xaemc", i, "_k_vec", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data_new = as.vector(data)
  assign(nam_new, data_new)
  temp_tdum_k_all[,j] = data_new
  rm(data)

  nam_precip = paste("precip_xaemc", i, "_kdat", sep = "")
  nam_new_precip = paste("precip_xaemc", i, "_k_vec", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_new_precip = as.vector(data_precip)
  assign(nam_new_precip, data_new_precip)
  precip_tdum_k_all[,j] = data_new_precip
  rm(data_precip)

  j = j + 1
}


for (i in Exp_list_D) {
  nam = paste("temp_xaemd", i, "_kdat", sep = "")
  nam_new = paste("temp_xaemd", i, "_k_vec", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data_new = as.vector(data)
  assign(nam_new, data_new)
  temp_tdum_k_all[,j] = data_new
  rm(data)

  nam_precip = paste("precip_xaemd", i, "_kdat", sep = "")
  nam_new_precip = paste("precip_xaemd", i, "_k_vec", sep = "")
  data_precip = matrix(data=eval(parse(text=nam_precip)), ncol = 201)
  data_new_precip = as.vector(data_precip)
  assign(nam_new_precip, data_new_precip)
  precip_tdum_k_all[,j] = data_new_precip
  rm(data_precip)

  j = j + 1
}

rm(list=ls()[grep("_vec",ls())])
rm(list=ls()[grep("kdat",ls())])

# --------------------------------------------------------------
# EXTRACT VARIABLES FOR EMULATOR (X AND Y)
# --------------------------------------------------------------

model_input_tdum = data.matrix(cont_paramdat)
model_output_temp = temp_tdum_k_all
model_output_precip = precip_tdum_k_all


# and put then on a grid format
model_output_temp <- array(model_output_temp, c(201, 201, 100 ))
model_output_precip <- array(model_output_precip, c(201, 201, 100 ))

rm(list= ls()[!(ls() %in% c('model_input_tdum','model_output_temp','model_output_precip'))])

# -----------------------------------------------------------------------------
# write all temperatures from the precursor climate model runs to NetCDF file
# -----------------------------------------------------------------------------
month="aug"
dims_lon = ncdim_def("lon","degrees_north",seq(1:201))
dims_lat = ncdim_def("lat","degrees_north",seq(1:201))
dims_exp = ncdim_def("exp_number","number",seq(1:100))
temp_all = ncvar_def("model_output_temp","degr K",list(dims_lon,dims_lat,dims_exp))
prefix="temp_"
suffix="_all_aism.nc"
ncnew    = nc_create(paste(prefix,month,suffix,sep=""),temp_all)
ncvar_put(ncnew,temp_all,model_output_temp)
nc_close(ncnew)

# -----------------------------------------------------------------------------
# write all precipitation files from the precursor climate model runs to NetCDF file
# -----------------------------------------------------------------------------
dims_lon = ncdim_def("lon","degrees_north",seq(1:201))
dims_lat = ncdim_def("lat","degrees_north",seq(1:201))
dims_exp = ncdim_def("exp_number","number",seq(1:100))
precip_all = ncvar_def("model_output_precip","degr K",list(dims_lon,dims_lat,dims_exp))
prefix="prec_"
suffix="_all_aism.nc"
ncnew    = nc_create(paste(prefix,month,suffix,sep=""),precip_all)
ncvar_put(ncnew,precip_all,model_output_precip)
nc_close(ncnew)
