source("pca_emul.R") # needs pca_reconstruct_mean

stochastic_pca_onestep  <-  function (PE, old.input, old.state, new.input) {
  # PE is the output of pe_c. This is the pca emulator
  # old.input is the input at time step n-1
  # new.input is the input at time step n
  # old.state is the output of this routine at time step n-1. 
  # old.input must be NULL if this is the first time step, in which case old.state is not read
  # this will generate one output made of a list composed of : simulated.field, which can
  # be passed on to the ice sheet model to compute the next ice state, and 
  # simulated.scores, which will be used in the next call to stochastic_pca_onestep
  
  if (is.null(old.input)) return( .stochastic_pca_onestep_null_oldinput ( PE, new.input ) )

  input_matrix  <-  matrix ( c(old.input, old.input), ncol = length(old.input), rnow = 2, byrow=TRUE)
  # input_matrix must have two rows, corresponding to the two vectors of input parameters
  GPp <- lapply(PE$GPList, function(EM_Cali)
                GP_P(EM_Cali, input_matrix, calc_var = TRUE, extra_output = FALSE ) )
  # each emulator generates a 'means' and 'variances' vector
  means     <- lapply(GPp, function(x) x$yp)
  variances <- lapply(GPp, function(x) x$Sp)
  simulated.pcscores <- mapply(
        function(yp, Sp, oldstate) {
        sx  <-  sqrt(Sp [1,1])
        sy <-   sqrt(Sp [2,2])
        rho  <- Sp[1,2]/sx/sy
      
        # mux and muy are the  _expected_ values of  scores at timestep i-1 and i
        mux <- yp[1]
        muy <- yp[2]
 
        mu_temp2  <-  muy + sy*rho * ( oldstate - mux ) / sx
        sd_temp2 <-  sy * sqrt ( 1 - rho*rho )
 
        return(rnorm ( mean = mu_temp2 sd = sd_temp2) )
        }, 
        # the following lists are the list of which the inputs
        # will be passed to the function which we have just defined
        means,variances, old.state$score_list)
  
   # now that we have simulated pc_scores we can reproduce the map
   
   simulated.field <- pca_reconstruct_mean (PE$L, simulated.pcscores)
   # and we still need to add a random number of residual  variance
   simulated.field  <-  simulated.field + 
                        rnorm(length(simulated.field), sd=sqrt(PE$L$rvm))
   
   return (list ( simulated.field = simulated.field, score_list =  simulated.pcscores )) }


 # this is a dummy function that the user should not call directly
# it is called by stochastic_pca_onestep when the old.input is null. 
  .stochastic_pca_onestep_null_oldinput <- function ( PE, new.input ) {
    GPp <- lapply( PE$GPList, function(EM_Cali)
                   GP_P(EM_Cali, new.input, calc_var = FALSE,
                        extra_output = FALSE ) )
    # each emulator generates a 'means' and 'variances' vector
    means     <- sapply(GPp, function(x) x$yp)
    variances <- sapply(GPp, function(x) x$Sp_diag)
    simulated.pcscores <- norm(length(means) , mean = means, 
                               sd = sqrt ( variances ) )
    simulated.field <- pca_reconstruct_mean (PE$L, simulated.pcscores)
    simulated.field  <-  simulated.field + rnorm(length(simulated.field), 
                                                 sd=sqrt(PE$L$rvm))
    return (list ( simulated.field = simulated.field, score_list =  simulated.pcscores )) 
  }
 
