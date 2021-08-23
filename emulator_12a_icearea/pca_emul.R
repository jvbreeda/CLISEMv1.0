# 
# Copyright (c) 2014 Michel Crucifix
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


#jvb 11022019: GP scripts added to this file
meanrm <- function (x) mean (x, na.rm=TRUE)

pca <-
function(A)
{
  nlon = dim(A)[1]
  nlat = dim(A)[2]
  nr = dim(A)[3]

  A = matrix(A, nlon*nlat, nr)

  means = apply(A,1,meanrm)
  Am = sweep( A,1, means)

  o = which( ! is.na (means ) )

  U = svd(Am[o,])

  O = A*NA
  O[o,] = U$u

  O = array(O, c(nlon,nlat,nr))

  # reformatting

  means = matrix(means, nlon, nlat)

  d     = U$d

  return(list(mean=means, PCA=O, amps = U$v, d=d))

}

funcmu_linear <-
function(x) 
				{x = c(1, x  )  
				return(x)
				}

funcmu_quadratic <-
function(x) 
				{ 
          out = funcmu_linear(x)
          n = length(x)
          for (i in seq(n)) for (j in seq(i,n)) 
          out = c(out, x[i] * x[j])
	      	       	    return(out)
						}

# overloads exp function for creation of sparse matrices
# commented out on
# Fri Oct  4 14:33:38 CEST 2013 as at the end
# we won't use sparse matrices
# exp <- function(x) Vectorize({ifelse(x<(-2), 0, .Primitive("exp")(x)) })

matern <- function(nu)
{
  require(gsl)
  if (nu > 100) return(exp)

  covar <- function(x)
  # uses will normall supply the minus squared distance
  {
    d = sqrt(-x)*2.
    if (d < 1.e-5) return (1)
    dn =  sqrt(nu)*d
    2^(1-nu)/ gsl_sf_gamma(nu) * dn^nu * bessel_Knu ( nu, dn )
  }

  Vectorize(covar)
}

cov_mat = function(lambda=lambda,X1=X1,X2=X2, covar=exp)
{
  # requires lambda to be a vector
  theta = lambda$theta
  nk = length(lambda$theta)
  nx = nrow(X1)
  ny = nrow(X2)
  RR = array(0, c(nx,ny,nk))
  for ( k in seq(nk) ) { RR[,,k] = outer(X1[,k], X2[,k], "-") / theta[k]  }
  R = covar(- ( apply(RR^2,c(1,2),sum) ))
}

cov_mat_1 = function(lambda=lambda,X1=X1,X2=X2, covar=exp)
{
  theta = lambda$theta
  # here need to adda check of consistency between theta and ncol X1 and nelem of X2
  nk = length(theta)
  nx = nrow(X1)
  vx = matrix(0,nx,nk)
  for  ( k in seq(1,nk) )
  {
    vx[,k] =  ((X1 [,k] - X2[k])/theta[k])^2
  }
  R = covar ( - ( apply(vx,1,sum)) )
}





GP_C <-
function( X, Y ,lambda, regress='linear', covar=exp )
  # revision history
  # 4.10.2013 : added covar option + passed in output.
  #             Backward campatible.
{
  if (is.character(regress))
  {
  funcmu = get(sprintf('funcmu_%s',regress))
  if  ( ! is.function(funcmu)) stop ('invalid regression model')
  if  ( ! is.matrix(Y)) Y=t(t(Y))
  if  ( ! is.matrix(X)) X=t(t(X))
  }
  else if (is.function(regress))
  {
  if (is.function(regress)) funcmu = regress
  }
  else stop('regress must either be a string or a function')

  muX  = t(apply(X, 1, funcmu))
  # note : if only a single row gets out of this this
  # probably means that muX is in fact a constant
  if ( nrow(muX) == 1) muX = (t(muX))

  nq   <- ncol(muX)

  n    <- nrow(X)
  nn   <- ncol(X)
  nbr  <- n - nq
  nbrr <- n - nq - 2

  R   <- cov_mat(lambda,X,X,covar)
  # R1X <- solve(R,muX) # for P matrix, disgarding the nugget

  # apply nugget (a la Andrianakis et Challenor)
  Rt   <- R + diag(n) * lambda$nugget
  R1tX <- solve(Rt,muX) # for P matrix

  dummy1 =  t(muX) %*% R1tX
  K      =  solve( dummy1, t(muX))
  betahat = K %*% solve(Rt, Y)


  dummy2  <- Y - muX %*% betahat
  e       <- solve(Rt, dummy2)
  # e is the notation of Oakley and OHagan 2004
  # equivalent to Nabila's formulation
  # See Bastos and O'Hagan 2009
  sigma_hat_2 <-  t( dummy2 ) %*%  e / nbrr

  # mean error at design points (proportionnal to nugget)

  M      <- (lambda$nugget)^2  * t(dummy2) %*% solve ( t(Rt) %*% Rt , dummy2 ) / n

  # max error at design points (when nugget tends to infinity, then just regression)

  Minfty <-  t(dummy2) %*% dummy2  / n
  # e.g. Baston and OHagan, eq. 15
  # Note : Adrianakis and Challenor have the same
  # result but multiplied by nbrr
  # (their definition of sigma_hat has this additional
  # factor, which propagates down to the likelihood)

  log_REML = -  1/2*( (nbr)*log(diag(sigma_hat_2))
                + log(det(Rt)) + log( det(t(muX) %*% R1tX )) )

  # penalised log lik (eq. 17 of Ardianakis and Challenor )
  # epsilon defaults to one


  if (is.null(lambda$epsilon)) lambda$epsilon = 1
  log_pen_REML = log_REML - 2 * (M / Minfty / lambda$epsilon )


  EM_Cali = list(betahat=betahat, sigma_hat_2=sigma_hat_2,
                 R=R, Rt = Rt,  muX = muX, X=X, Y=Y, lambda=lambda, e=e,
                 funcmu=funcmu, # R1X = R1X,
                 R1tX=R1tX , log_REML = log_REML,
                 log_pen_REML=log_pen_REML, covar = covar, nbrr=nbrr )

  attr(EM_Cali, "class") <-  "GP_Emul"
  return(EM_Cali)
}


logLik.GP_Emul <- function(E)
  {
    return(E$log_REML)
  }


BIC.GP_Emul <- function(E)
  {
    k <- with(E, length(lambda$theta) + length(lambda$nugget))
    n <- with(E, length(Y) )
    bic <- with (E,  -2 * log_REML +  k * (log ( n ) - log ( 2 * pi ) )  )
    return(bic)
  }



GP_P <-
function( EM_Cali, x, calc_var=FALSE, extra_output=FALSE) 
{
  # revision history
  # 4.10.03 : added covar
  if (! is.matrix(x)) x=as.matrix(x)
  X <- EM_Cali$X
  Y <- EM_Cali$Y
  # the Rt includes the nugget 
  R <- EM_Cali$R
  e <- EM_Cali$e
  Rt <- EM_Cali$Rt
  #R1X  <- EM_Cali$R1X
  R1tX <- EM_Cali$R1tX
  covar  <- EM_Cali$covar

  lambda <- EM_Cali$lambda
  betahat   <- EM_Cali$betahat
  sigma_hat_2 <- EM_Cali$sigma_hat_2
  funcmu       <- EM_Cali$funcmu
  muX  <- EM_Cali$muX
  mux  <- t(apply(x, 1, funcmu)) 
  
  # if only a constant in returned need to transpose output
  if (length(funcmu(1)) == 1) {mux = t(mux)}

  nx   <- nrow(x)
  n    <- nrow(X)    
  # check if Nabila's code consistent with this. Was previousy
  # (in error) ncol(X)
  nq   <- length(betahat)
  nbrr <- n - nq - 2 

  # r is T(x) in the challenor paper. No nugget ! 
  r  <- if (nx == 1) { matrix(cov_mat_1(lambda,X,x,covar),n,1)
                     } else {cov_mat(lambda,X,x,covar)}

  yp_mean = mux %*% betahat 
  yp_gaus = t(r) %*% e

  yp = yp_mean + yp_gaus 

  if (calc_var) # do we need the full output covariance matrix ? 
  {
    rr <-  if ((nx)==1) {1}   else {cov_mat(lambda,x ,x,covar)}

    # save a copy of nugget-free covariance
    rr_nuggetfree <- rr

    rr <- rr + diag(nx) * lambda$nugget
    # correction : the nugget must not only be added to the diag, but in fact
    # to any couple of inputs that would be identical 

    Ind = cbind(which(duplicated(x)),nrow(x)+1-which(duplicated(apply(x, 2, rev) )))
    Ind = rbind(Ind, Ind[, c(2,1)])
    rr[Ind] = rr[Ind] + lambda$nugget 



    # see email to Adrianikis and Challenor 5.08.2013 without response
    # whether one hould use R1X (i. e. : A-1) or R1tx (i. e. Atilde - 1)
    # new email on 22 august
    # ok. response on 9 septembre : Must use Atilde. 
    # equation modified accordingly
    P  <- ( mux - t(r) %*% R1tX )

    ## in the notation Oakley OHagan : 

   
    cxx0 <- -  t(r) %*% solve ( Rt, r ) + ( P  %*% solve ( (t(muX) %*% R1tX )  ,  t(P) ) )

    cxx_star = rr  + cxx0
    cxx_star_nuggetfree = rr_nuggetfree  + cxx0

    Sp = cxx_star * as.numeric(sigma_hat_2)
    Sp_nuggetfree = cxx_star_nuggetfree * as.numeric(sigma_hat_2)

    Sp_diag = diag(Sp)
    Sp_diag_nuggetfree = diag(Sp_nuggetfree)

    if (extra_output)
    {
      # extra output useful to compute the triple integrals fro
      # general sensitivity analysis after Oakley and Ohagan
      # attention : our ht is the transpose of their 'h', but
      # easier for the calculations that will follow

      # output sparse matrices for cxx, htt, ttt
      # to do : this needs to be an option
      # i believe that for hht this may be counterproductive
      # or also if lambda is large so that in practice there are
      # little excluded points, or at last if the 
      # user does not cut exponential tails
      # to do: select sparse matrix or dense matrix
      # depending on th proportion of zeros in r. 
      
      # extra note on Fri Oct  4 14:04:32 CEST 2013
      # it seems that all of this is actually a waste of 
      # computing time : theses output are not used for
      # varanal_fast (which was verified to produce the same
      # output as varanal). 
 
      mux = as.matrix(mux)
      r = as.matrix(r)
      Emul_pred = list(yp=yp, yp_mean=yp_mean, yp_gaus=yp_gaus, Sp=Sp, Sp_nuggetfree = Sp_nuggetfree, Sp_diag = Sp_diag, Sp_diag_nuggetfree = Sp_diag_nuggetfree, r=r, ht=mux, 
                       cxx = rr, 
                       cxx_star = cxx_star,
                       hht = aperm ( outer(t(mux), mux) , c(2,3,1,4) ),
                       htt = aperm(outer(t(mux), t(r) ), c(2,3,1,4)) , 
                       ttt = aperm(outer(r, t(r) )  , c(2,3,1,4))   
                       )
    }
    else 
    Emul_pred = list(yp=yp, Sp=Sp, Sp_diag = Sp_diag, 
     Sp_diag_nuggetfree, yp_mean=yp_mean, yp_gaus=yp_gaus)

  } else
  {
    rr_nuggetfree = 1.
    rr =  rr_nuggetfree + lambda$nugget
    dummy1 = (t(muX) %*% R1tX ) 
    Sp_diag=rep(0,nx)
    Sp_diag_nuggetfree=rep(0,nx)
    for (j in seq(1,nx))
      {
        rj = r[,j, drop=FALSE]
        # bug corrected on Fri Nov 14 22:48:46 CET 2014
        # R1X below -> R1tX
        P  <- ( mux[j,] - t(rj) %*% R1tX )
        cxx0 = - (t(rj) %*% solve(Rt, rj)) + P %*%  solve ( dummy1 , t(P) )  
        cxx = rr + cxx0
        cxx_nuggetfree = rr_nuggetfree + cxx0
        Sp_diag[j] = sigma_hat_2 * cxx 
        Sp_diag_nuggetfree[j] = sigma_hat_2 * cxx_nuggetfree
      }
    if (extra_output) 
      Emul_pred = list(yp=yp, Sp_diag=Sp_diag, Sp_diag_nuggetfree = Sp_diag_nuggetfree, ht=mux, tt=t(r))
    else 
      Emul_pred = list(yp=yp, Sp_diag = Sp_diag,  Sp_diag_nuggetfree = Sp_diag_nuggetfree, yp_mean=yp_mean, yp_gaus=yp_gaus)
  }

}


#end of GP scripts



mypca <- function(var, nkeep)
{
    pca = pca(var)
  # rvm : residual variance estimated  at each grid point
  pca$rvm = sum ( pca$d[-seq(nkeep)]^2) / length(which( ! is.na (var[,,1]))) / dim(var)[3]
    pca
}


pca_reconstruct_mean <- function(L,means)
{
  # sum over basis elements scaled by L$d + mean
  sv = seq(length(means))
  nv = length(means)
  return ( apply ( sweep(L$PCA[,,sv], 3, means * L$d[sv] ,'*' ) , c(1,2), sum)  + 
     L$mean )
}

pca_sample_quantiles <- function (L, means, vars, dof, n = 10000)
{
  nv = (length(means))
  sv = seq(length(means))
  # sample variances
  #m = outer ( sqrt ( vars ) , rt(n , dof, 0)   ) 
  m = matrix( rt(n * nv , dof, 0) , nv, n)
  m = sweep ( m, 1, sqrt ( vars) , '*') 
  m = sweep ( m, 1, means, '+') 
  m = sweep ( m, 1, L$d[sv], '*' ) 

 # add means

  # provisionnaly reshapes LL to prepare the inner product
  no = dim(L$PCA)[1] * dim(L$PCA)[2]
  LL = matrix ( L$PCA[,,sv], no, nv)
  # multiply by means
  R  = LL %*% m 

  # residal variance
  rv = rnorm(n , 0, sd =  sqrt ( L$rvm ) )
  R = R + rv 

  # reshape under friendly form
  R  = array ( R, c(dim(L$PCA)[1], dim(L$PCA)[2],  n ) ) 
  R  = sweep(R, c(1,2), L$mean, '+')

 

  # quantiles                                           1    2      3    4     5     6      7    Inf
  O = apply( R, c(1,2), function (i) quantile(i, probs=c(0.002, 0.025, 0.17, 0.5, 0.83, 0.975, 0.998) ) ) 
  return(O)
}

compare_quantiles_with_ref <- function ( O, REF)
{
  I = sweep(O,  c(2,3) , REF) 
  Q = apply(I, c(2,3), function (i) min ( which ( i > 0 )) ) 
  q1 = length(which(Q == 4 | Q == 5 ))
  q2 = length(which(Q == 3 | Q == 6 ))
  q3 = length(which(Q == 2 | Q == 7 ))
  q4 = length(Q) - q1 - q2 - q3
  DF = sqrt ( mean  (REF - O[4,,] , na.rm = TRUE)^2 ) 
  return(c(q1,q2,q3,q4, DF))
}


pca_reconstruct_var <- function(L,variances)
{
  # sum of component variances provided by variances 
  # n is number of experiments used to estimate covariance
  sv = seq(length(variances))
  nv = length(variances)
  var=  apply ( sweep(L$PCA[,,sv]^2, 3, variances * (L$d[sv]^2) ,'*' ) , c(1,2), sum)  
  # spread residual variance equally
  # may not be a good idea: better to spread variance according 
  # to residual space
  if (! is.null (L$rvm)) 
  {
  var  <- var + L$rvm
  }
  return ( var ) 
}

pe_l1o <- function (X, Y, pca_function=pca, hp, nkeep=10, ...) 
{
  # X : input vector
  # Y : outputs (as matrix)
  # 

  VAR= lapply( seq(nrow(X)) , function (i) {
    PE  = pe_c (X[-i,], Y[,,-i], pca_function=pca_function, hp=hp,nkeep=nkeep)
    PRE = pe_p (X[i,, drop=FALSE], PE)
    # the distance metric is here blindly estimated assuming independent pixels
    DIFF   =  - PRE$mean + Y[,,i]
    DS  =  abs (PRE$mean - Y[,,i] ) / sqrt (  PRE$var )
    list(DIFF=DIFF, DS=DS) })
    VAR
}


pe_l1o_barplot_student <- function (X, Y, pca_function=pca, hp, nkeep=20, ...) 
{
  # X : input vector
  # Y : outputs (as matrix)
  # 

  VAR= sapply( seq(nrow(X)) , function (i) {
    PE  = pe_c (X[-i,], Y[,,-i], pca_function=pca_function, hp=hp,nkeep=nkeep)
    OUT = pe_p_student(X[i,, drop=FALSE], PE)

    compare_quantiles_with_ref(OUT, Y[,,i] )
    })

    # tyding up : set DF as an attribute
    DF = VAR[5,]
    VAR = VAR[1:4,]
    attr(VAR, 'DF') <- DF
    VAR
  
}

pe_l1o_maxerr <- function (X, Y, pca_function=pca, hp, nkeep=20, ...) 
{
 VAR= lapply( seq(nrow(X)) , function (i) {
    PE  = pe_c (X[-i,], Y[,,-i], pca_function=pca_function, hp=hp,nkeep=nkeep)
    PRE = pe_p (X[i,, drop=FALSE], PE)
    # the distance metric is here blindly estimated assuming independent pixels
    DS  =  abs (PRE$mean - Y[,,i] ) / sqrt (  PRE$var )
    })
    print('ici')
    VAR = simplify2array(VAR)
    apply(VAR, c(1,2), max)
}


pe_l1o_barplot <- function (X, Y, pca_function=pca, hp, nkeep=20, ...) 
{
  # X : input vector
  # Y : outputs (as matrix)
  # 

  VAR= sapply( seq(nrow(X)) , function (i) {
    PE  = pe_c (X[-i,], Y[,,-i], pca_function=pca_function, hp=hp,nkeep=nkeep)
    PRE = pe_p (X[i,, drop=FALSE], PE)
    # the distance metric is here blindly estimated assuming independent pixels
    DF   = sqrt ( mean  (PRE$mean - Y[,,i] , na.rm = TRUE)^2 ) 
    DS  =  abs (PRE$mean - Y[,,i] ) / sqrt (  PRE$var )
    q1 <- length ( which(DS < 1))
    q2 <- length ( which(DS < 2)) - q1
    q3 <- length ( which(DS < 3)) - q2 - q1
    q4 <- length ( which(DS < 99999.)) - q3 - q2 - q1 
    c(q1,q2,q3,q4,DF)})

    # tyding up : set DF as an attribute
    DF = VAR[5,]
    VAR = VAR[1:4,]
    attr(VAR, 'DF') <- DF
    VAR
  
}

#mypca was pca jvb 11022019
#pe_c <- function (X, Y, pca_function=mypca, hp, nkeep)
pe_c <- function (X, Y, pca_function=mypca, hp, nkeep, ...)
{
  # emulator calibration (for same l and nugget applied to all pcas ! )
  # otherwise one may look at an autmatic calibration (later)
  # X : design
  # Y : output under the form (nlat,nlon, nr) where nr is the
  #     number of design members
  # hp : table of hyperparameters
  # output : generate pca and list of gaussian process emulators using
  # pca funcion and table of hyperparameters

  L <- pca_function(Y, nkeep)
  n <- nrow(X)
  m <- ncol(X)

  L$scaled_amps <- sweep( L$amps, 2  , L$d, '*') 
  # attention : temporary fix : set residual variance to zero !!! 
  # do not forget to change
  GPx <- lapply(seq(nkeep), function(i) 
              GP_C(X[,], L$amps[,i], 
                   list(theta=hp[1:m,i], nugget=hp[m+1,i])) )

  return(list(L=L, GPList = GPx))

}

pe_p <- function (x, PE)
{
  # pca emulator predictor, at input x
  GPp <- lapply(PE$GPList, function(EM_Cali) GP_P(EM_Cali, x, calc_var = FALSE, extra_output = FALSE ) )
  means     <- sapply(GPp, function(x) x$yp)
  variances <- sapply(GPp, function(x) x$Sp_diag)  
  # need to consider residual variance !!! 
  # number of samples used for covariance estimation
  n = nrow(PE$GPList[[1]])
  mean_field <- pca_reconstruct_mean ( PE$L, means ) 
  var_field  <- pca_reconstruct_var ( PE$L, variances)
  return(list(mean=mean_field, var=var_field, means=means, variances=variances))
}

pe_p_student <- function (x, PE)
{
  # pca emulator predictor, at output x
  GPp <- lapply(PE$GPList, function(EM_Cali) GP_P(EM_Cali, x, calc_var = FALSE, extra_output = FALSE ) )
  # should figure out dof here
  dof = 20
  means     <- sapply(GPp, function(x) x$yp)
  variances <- sapply(GPp, function(x) x$Sp_diag)  
  out <- pca_sample_quantiles(PE$L, means, variances, dof)
  # need to consider residual variance !!! 
  return(out)
}





