# --------------------------------------------------
# LOAD LIBRARIES
#library(ncdf)
#jvb
library(ncdf4)
library(fields)
library(maps)
#library(mapproj)
# SET WORKING DIRECTORY

# to set a working directory
#setwd("/Users/jvbreeda/Desktop/PhD/HadCM3/matlab/hadcm3_output/eocene/xaenaa/pd/40km")

#Exp = "tdum"
Exp1="xaema"
Exp2="xaemb"
Exp3="xaemc"
Exp4="xaemd"
month="jan"
prefix="mean_pd_"
suffix=".nc"

Exp_list_A = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
#
Exp_list_B = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")

Exp_list_C = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
#
Exp_list_D = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v")
# -------------------------------------------------
# READ IN INPUT DATA
cont_paramdat = read.table('emul_input_5variables.txt', header=TRUE)
cont_param_dim=dim(cont_paramdat)
obl = cont_paramdat$obliquity
esinw = cont_paramdat$esinw
ecosw = cont_paramdat$ecosw
co2 = cont_paramdat$co2_ppm
#ecc = cont_paramdat$eccentricity
#omega = cont_paramdat$omega
ice = cont_paramdat$ice


# READ IN OUTPUT DATA
temp=list((0))

    togrid=file.path('..','output_emulator_100b',"AISM40_eot.nc")
    ncgrid = nc_open(togrid,write=TRUE)
    lats_dat = ncvar_get(ncgrid, "lat")
    lons_dat = ncvar_get(ncgrid, "lon")


for (i in Exp_list_A) {
    tofile=file.path('..','output_emulator_100b',paste(Exp1,i,sep=""),paste(prefix,month,suffix,sep=""))
    nc = nc_open(tofile,write=TRUE)

	temp = ncvar_get(nc, "temp")
	nam = paste("temp_xaema", i, "_kdat", sep = "")
	assign(nam, temp)
	nc_close(nc)
	rm(nc, temp)
}

for (i in Exp_list_B) {
    tofile=file.path('..','output_emulator_100b',paste(Exp2,i,sep=""),paste(prefix,month,suffix,sep=""))

    nc = nc_open(tofile,write=TRUE)

      temp = ncvar_get(nc, "temp")
      nam = paste("temp_xaemb", i, "_kdat", sep = "")
      assign(nam, temp)
      nc_close(nc)
      rm(nc, temp)
}

for (i in Exp_list_C) {
    tofile=file.path('..','output_emulator_100b',paste(Exp3,i,sep=""),paste(prefix,month,suffix,sep=""))
    nc = nc_open(tofile,write=TRUE)

      temp = ncvar_get(nc, "temp")
      nam = paste("temp_xaemc", i, "_kdat", sep = "")
      assign(nam, temp)
      nc_close(nc)
      rm(nc, temp)
}


for (i in Exp_list_D) {
    tofile=file.path('..','output_emulator_100b',paste(Exp4,i,sep=""),paste(prefix,month,suffix,sep=""))
    nc = nc_open(tofile,write=TRUE)

      temp = ncvar_get(nc, "temp")
      nam = paste("temp_xaemd", i, "_kdat", sep = "")
      assign(nam, temp)
      nc_close(nc)
      rm(nc, temp)
}


# MAKE MISSING VALUES

for (i in Exp_list_A) {
  nam = paste("temp_xaema", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)
}


for (i in Exp_list_B) {
  nam = paste("temp_xaemb", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)
}


for (i in Exp_list_C) {
  nam = paste("temp_xaemc", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)
}


for (i in Exp_list_D) {
  nam = paste("temp_xaemd", i, "_kdat", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data = replace(data, data == 2.00000004008175e+20, NA)
  assign(nam, data)
  rm(data)
}


# --------------------------------------------------
# VECTORIZE TEMPERATURE MATRICES (BY COLUMN (LATITUDE))

temp_tdum_k_all = matrix(0, nrow = length(temp_xaemaa_kdat), ncol = 100)
j = 1


for (i in Exp_list_A) {
  nam = paste("temp_xaema", i, "_kdat", sep = "")
  nam_new = paste("temp_xaema", i, "_k_vec", sep = "")
  data = matrix(data=eval(parse(text=nam)), ncol = 201)
  data_new = as.vector(data)
  assign(nam_new, data_new)
  temp_tdum_k_all[,j] = data_new
  rm(data)
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
  j = j + 1
}


rm(list=ls()[grep("_vec",ls())])
rm(list=ls()[grep("kdat",ls())])

# EXTRACT VARIABLES FOR EMULATOR (X AND Y)

model_input_tdum = data.matrix(cont_paramdat)

model_output_tdum = temp_tdum_k_all

# and put then on a grid format

model_output_tdum <- array(model_output_tdum, c(201, 201, 100 ))

rm(list= ls()[!(ls() %in% c('model_input_tdum','model_output_tdum'))])

##################################
# write temp to NetCDF
##################################

dims_lon = ncdim_def("lon","degrees_north",seq(1:201))
dims_lat = ncdim_def("lat","degrees_north",seq(1:201))
dims_exp = ncdim_def("exp_number","number",seq(1:100))
temp_all = ncvar_def("model_output_tdum","degr K",list(dims_lon,dims_lat,dims_exp))
ncnew    = nc_create("temp_jan_all_aism.nc",temp_all)
ncvar_put(ncnew,temp_all,model_output_tdum)
nc_close(ncnew)
