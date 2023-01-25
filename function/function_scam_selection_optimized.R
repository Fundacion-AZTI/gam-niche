# Script Info -------------------------------------------------------------

# Function to perform forward model selection based on AIC and p-value for scam models with concave terms (IM-20-MISSION)

# The function arguments are: 
#   basef   : base model formula, usually  y ~ 1
#   vars    : list of potential variables to be included
#   dat     : data set
#   vmax    : maximum number of variables to be included. If NULL, there is no upper limit
#   aic.tol : stopping criterion. When the aic improvement is less than aic.tol, we stop the model selection process
#   k       : knots
#   sp      : either the fixed smoothing parameter or NULL (if sp are estimated) 
#   plevel       : if p=NULL, only AIC criteria will be used, if plevel is a number, the function will only include smooth whose p-value is smaller than plevel 

# Leire Ibaibarriaga, Mireia Valle and Guillem Chust (based on code from Leire Citores for IM-18-ANICHO)
# Last modifications 7/04/2022


# Function for scam models ------------------------------------------------

modsel.scam <- function(basef, vars, dat, vmax=NULL, aic.tol=2, k=8, sp= NULL, plevel=0.05){
  require("scam")
  # initialise the vector of selected variables
  svars <- NULL 
  # initialise the list of selected models at each step
  smod <- list() 
  # base model 
  m0 <- scam(formula(basef), family=binomial(link="logit"), data=dat)
  aic0 <- AIC(m0)
  smod[[1]] <- m0
  # maximum number of variables to be included
  vmax <- min(vmax, length(vars)) 
  # loop: add one variable at each step 
  for (i in 1:vmax){ 
    print(paste("Fitting models with", i, "terms"))
    m.list <- list() # list of models with i terms
    for (j in 1:length(vars)){
      if (! vars[j] %in% svars){
        tempvars <- c(svars, vars[j])
        formu <- paste(c(basef, paste0("s(",tempvars,",m=2,bs='cv',k=",k,")")), collapse="+")
        cat("\t", formu, "\n")
        if (is.null(sp)){
          m <- scam(formula(formu), family=binomial(link="logit"), data=dat, sp=NULL)
        }else{
          m <- scam(formula(formu), family=binomial(link="logit"), data=dat, sp=rep(sp,length(tempvars)))
        }
        m.list[[vars[j]]] <- m
      }
    }
    if (!is.null(plevel)){
      sel <- sapply(m.list, function(x) all(summary(x)$s.pv < plevel & !is.na(summary(x)$s.pv)))  # check that all the terms are significant
      if (sum(sel, na.rm=T)>0){
        m.list <- m.list[sel]
      }else{
        print(paste("The smooths were not significant for ", plevel))
        break
      }
    }
    aic.vec <- sapply(m.list, AIC)
    aic.min <- aic.vec[which.min(aic.vec)]
    if ( (aic0-unname(aic.min)) > aic.tol){
      svars <- c(svars, names(aic.min)) # selected vars with smallest AIC
      smod[[i+1]] <- m.list[[names(aic.min)]] # selected model
      aic0 <- unname(aic.min)
    }else{
      print(paste("The AIC did not improve in more than", aic.tol, "units"))
      break
    }
  }
  names(smod) <- c("null", svars)
  out <- list(smod=smod, svars=svars)
  return(out)
}


# End of script -----------------------------------------------------------


