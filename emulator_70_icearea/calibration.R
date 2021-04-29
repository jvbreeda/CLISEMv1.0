# main routine

# load resources
require(GP)
source('pca_emul.R')
require(BBmisc)


#------------------------------------------------------------------------------------
#		calibrate temperature
#------------------------------------------------------------------------------------

# calibrate the emulator as a function of the month
month=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

calibrate_emulator_temp <- function (month) {

# load the temperature fields from the climate model runs
prefix="xaen_aism_"
suffix=".R"
source(paste(prefix,month,suffix,sep=""))

X = model_input_tdum			 #model input: ecc, esinw, ecosw, co2, ice,...
Y = model_output_temp 			 #temperature
sx1=sum(X[,1])				 #esinw	
sx2=sum(X[,2])				 #ecosw
sx3=sum(X[,3])				 #obliquity
sx4=sum(X[,4])				 #CO2 concentration
sx5=sum(X[,5])				 #ice volume

#normalize the values before using them
b=X[,1]/sx1
c=X[,2]/sx2
d=X[,3]/sx3
e=X[,4]/sx4
f=X[,5]/sx5

#save the normalized output to use later for the prediction
save(sx1,file="sx1.Rda")
save(sx2,file="sx2.Rda")
save(sx3,file="sx3.Rda")
save(sx4,file="sx4.Rda")
save(sx5,file="sx5.Rda")

# make matrix from the normalized values
X = cbind(b[1:70], c[1:70], d[1:70], e[1:70], f[1:70]) 


# number of Principal Components to take into account
  nkeep=20 
# length scale and nugget hyperparameter values 
  hp = data.frame( l.esinw=rep(3.65, nkeep), l.ecosw = 3.31, l.obl = 5.22, l.co2 = 6, l.icevol = 1.2 , nugget = 0.0001)

# The way the pe_c routine takes hp is not okay: it requires a matrix shaped this way
  hp <- t( as.matrix(hp))

#load the emulator source code
source('./pca_emul.R')
E = pe_c (X, Y, mypca, hp=hp, nkeep=nkeep)
}

# calculate and save the mean and PCA scores for each month
E_list_temp <- lapply ( month, calibrate_emulator_temp)
save(E_list_temp,file="E_list_temp.Rda")


#------------------------------------------------------------------------------------
#		calibrate precipitation
#------------------------------------------------------------------------------------


month=c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

calibrate_emulator_precip <- function (month) {

# load the precipitation fields from the climate model runs
prefix="xaen_aism_"
suffix=".R"
source(paste(prefix,month,suffix,sep=""))

X = model_input_tdum			 #model input: ecc, esinw, ecosw, co2, ice,...
Y = model_output_precip			 #precipitation
sx1=sum(X[,1])				 #esinw	
sx2=sum(X[,2])				 #ecosw
sx3=sum(X[,3])				 #obliquity
sx4=sum(X[,4])				 #CO2 concentration
sx5=sum(X[,5])				 #ice volume

#normalize the values before using them
b=X[,1]/sx1
c=X[,2]/sx2
d=X[,3]/sx3
e=X[,4]/sx4
f=X[,5]/sx5

#save the normalized output to use later for the prediction
save(sx1,file="sx1.Rda")
save(sx2,file="sx2.Rda")
save(sx3,file="sx3.Rda")
save(sx4,file="sx4.Rda")
save(sx5,file="sx5.Rda")

# make matrix from the normalized values
X = cbind(b[1:70], c[1:70], d[1:70], e[1:70], f[1:70]) 

# calibrate the emulator as a function of the month
# number of Principal Components to take into account
nkeep=17
# length scale and nugget hyperparameter values
hp = data.frame( l.esinw=rep(0.5, nkeep), l.ecosw = 0.5, l.obl = 0.5, l.co2 = 0.5, l.icevol = 1.2 , nugget = 0.01)
hp <- t( as.matrix(hp))

#load the emulator source code
source('./pca_emul.R')
E = pe_c (X, Y, mypca, hp=hp, nkeep=nkeep)
}

# calculate and save the mean and PCA scores for each month
E_list_precip <- lapply ( month, calibrate_emulator_precip)
save(E_list_precip,file="E_list_precip.Rda")

