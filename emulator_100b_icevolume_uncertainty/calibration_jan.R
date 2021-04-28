# main routine
# calibrate the January Temperature

# load resources
require(GP)
source('pca_emul.R')
require(BBmisc)

prefix="xaem_aism_"
suffix=".R"
month="jan"


# load the January temperatures from the climate model runs
prefix="xaem_aism_"
suffix=".R"
source(paste(prefix,month,suffix,sep=""))

X = model_input_tdum		#model input: ecc, esinw, ecosw, co2, ice,...
Y = model_output_tdum		#temperature/precipitation

sx1=sum(X[,1])			#esinw	
sx2=sum(X[,2])			#ecosw
sx3=sum(X[,3])			#obliquity
sx4=sum(X[,4])			#CO2 concentration
sx5=sum(X[,5])			#ice volume

#normalize the values before using them
b=X[,1]/sx1
c=X[,2]/sx2
d=X[,3]/sx3
e=X[,4]/sx4
f=X[,5]/sx5


save(sx1,file="sx1.Rda")
save(sx2,file="sx2.Rda")
save(sx3,file="sx3.Rda")
save(sx4,file="sx4.Rda")
save(sx5,file="sx5.Rda")

# make matrix from the normalized values
X = cbind(b[1:100], c[1:100], d[1:100], e[1:100], f[1:100]) 

# number of Principal Components to take into account
nkeep=20

# length scale and nugget hyperparameter values 
hp = data.frame( l.esinw=rep(3.65, nkeep), l.ecosw = 3.31, l.obl = 5.22, l.co2 = 3.5, l.icevol = 1.2 , nugget = 0.001)

# The way the pe_c routine takes hp is not okay: it requires a matrix shaped this way
hp <- t( as.matrix(hp))

#load the emulator source code
source('./pca_emul.R')

# calibrate the GP PCA emulator with the routine pe_c
E_list_jan2 = pe_c (X, Y, mypca, hp=hp, nkeep=nkeep)

#save E_list and call it again
save(E_list_jan2,file="E_stochastic.Rda")

