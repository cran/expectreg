laws <-
function(B,DD,yy,pp,lambda,smooth,nb,nbp,nbunp,center,types,LAWSmaxCores=1)
    ###### expectile regression according to eilers, schnabel
    # parameters:
    # formula - vector of responses ~ f(vector of independent,type="which base to use") + ...
    # smooth - if smoothing with schall's algorithm, asymmetric cross validation or no smoothing shall be done
    # lambda - smoothing penalty, important if no smoothing is done
    # nb - number of parameters per covariate
    # nbp - number of penalized parameters per covariate
    # nbunp - number of unpenalized parameters per covariate
{
    nterms = length(nb)
    m = length(yy)
    np = length(pp)
    
    
    #  myapply <- lapply
    #  if (.Platform$OS.type == "unix" && require("parallel")) 
    #  {
    #      if (!parallel:::isChild()) 
    #      {
    #          myapply <- mclapply
    #      }
    #  }
    
    if(length(lambda) < nterms)
        lala = rep(lambda[1],nterms)
    else
        lala = lambda
    
    lala_orig <- lala
    
    dummy.reg <- function(pp,lala,smooth,yy,B,DD,nb,nterms,center)
    {
        #cat("Expectile: ",pp,"\n")
        
        if(smooth == "schall")
        {
            glatterms <- which(unlist(types) != "parametric" )
            # How many parameters are used after each covariate?
            nbc <- cumsum(c(1,nb)) # c(1,nb): Number of paramters for each covariate incl. intercept
                                   # cumsum(...): How many parameters have been used for each parameter at it's end
            
            # Shift of meaning: How many parameters are used before each covariate: Intercept vanishes and last term is irrelevant
            nbpc <- nbc + 1 + c(nbunp,0) # nbc + 1: Position of starting point for next covariate
                                         # ... + nbunp: shift of starting point (schall) due to unpenalized parts in this covariate.
            
            #cat(glatterms,"\n",nb,"\n",nbunp,"\n",nbp,"\n",nbc,"\n",nbpc,"\n\n\n\n\n\n\n\n\n\n")
            
            sch = schallCPPfun(glatterms = glatterms, y = as.numeric(yy), B = B, 
                               tau = as.numeric(pp), lambdashort_in = as.numeric(lala), 
                               DD = DD, NB=as.integer(nb), NBP=as.integer(nbp), 
                               NBPC=as.integer(nbpc), center=center)
            
            #cat(glatterms,"\n",nb,"\n",nbunp,"\n",nbp,"\n",nbc,"\n",nbpc,"\n\n\n\n\n\n\n\n\n\n")
            
            lala = sch$lambdashort
            vector.a.ma.schall = sch$a
            diag.hat = sch$diag.hat.ma
            if(sch$it > 99) {
                warning("Schall algorithm did not converge. Stopping after 100 iterations.")
            }
            if(sch$iter > 49) {
                warning("IWLS weights did not converge after 50 iterations.")
            }
        }
        else if(smooth == "aic"||smooth == "bic"||smooth == "ocv"||smooth == "gcv")
        {
            lala_glatt <- lala[unlist(types) != "parametric"]
            glatterms <- which(unlist(types) != "parametric" )
            
            if(length(glatterms) > 0) {
                acv.min = nlminb(start=lala_glatt,objective=smoothCPPP,y=as.numeric(yy),B=B,tau=pp,lambdashort_orig=as.numeric(lala_orig),DD=DD,NB=as.integer(nb),glatterms = glatterms, smoothtype=smooth,lower=0,upper=10000)
                lala[unlist(types) != "parametric"] <- abs(acv.min$par)
            }

            aa <- lltLS(y=as.numeric(yy), B=B, tau=pp, lambdashort_glatt=as.numeric(lala), lambdashort_orig=as.numeric(lala), DD=DD, NB=as.integer(nb), glatterms = 1:length(lala))
                            
            vector.a.ma.schall <- aa$a  
            diag.hat = aa$diag.hat.ma
        }
        
        else if(smooth == "cvgrid")
        {
            lala_glatt <- lala[unlist(types) != "parametric"]
            glatterms <- which(unlist(types) != "parametric" )
            
            
            if(length(glatterms) > 0) {
            lala_glatt = cvgrid_fast(yy,B,pp,DD,as.integer(nb),types,lala)
            lala[unlist(types) != "parametric"] <- abs(lala_glatt)
            }

            aa <- lltLS(y=as.numeric(yy), B=B, tau=pp, lambdashort_glatt=as.numeric(lala), lambdashort_orig=as.numeric(lala), DD=DD, NB=as.integer(nb), glatterms = 1:length(lala))
                            
            vector.a.ma.schall <- aa$a  
            diag.hat = aa$diag.hat.ma
        }
        else if(smooth == "lcurve")
        {
            lala_glatt <- lala[unlist(types) != "parametric"]
            glatterms <- which(unlist(types) != "parametric" )
            
            if(length(glatterms) > 0) {
            lala_glatt = lcurve_fast(yy,B,pp,DD,nb,types,lala)
            lala[unlist(types) != "parametric"] <- abs(lala_glatt)
            }

            aa <- lltLS(y=as.numeric(yy), B=B, tau=pp, lambdashort_glatt=as.numeric(lala), lambdashort_orig=as.numeric(lala), DD=DD, NB=as.integer(nb), glatterms = 1:length(lala))
                            
            vector.a.ma.schall <- aa$a  
            diag.hat = aa$diag.hat.ma
            
        }
        else
        {
            lala_glatt <- lala[unlist(types) != "parametric"]
            glatterms <- which(unlist(types) != "parametric" )
            
            aa <- lltLS(y=as.numeric(yy), B=B, tau=pp, lambdashort_glatt=as.numeric(lala), lambdashort_orig=as.numeric(lala), DD=DD, NB=as.integer(nb), glatterms = 1:length(lala))
            
            vector.a.ma.schall <- aa$a  
            diag.hat = aa$diag.hat.ma 
            
        }
        
        list(vector.a.ma.schall,lala,diag.hat)
    }
    
    if (.Platform$OS.type == "unix")
        coef.vector = mclapply(pp,function(pp) dummy.reg(pp,lala,smooth,yy,B,DD,nb,nterms,center),mc.cores = max(1,min(detectCores()-1,LAWSmaxCores)))
    else if (.Platform$OS.type == "windows")
        coef.vector = mclapply(pp,function(pp) dummy.reg(pp,lala,smooth,yy,B,DD,nb,nterms,center),mc.cores = 1)
    
    
    lala <- matrix(lambda, nrow=nterms, ncol=np)
    vector.a.ma.schall <- matrix(NA, nrow=sum(nb)+(1*center),ncol=np)
    diag.hat = matrix(NA,nrow=m,ncol=np)
    
    
    for(i in 1:np)
    {
        vector.a.ma.schall[,i] = coef.vector[[i]][[1]]
        lala[,i] = coef.vector[[i]][[2]]
        diag.hat[,i] = coef.vector[[i]][[3]]
    }
    
    return(list(vector.a.ma.schall,lala,diag.hat))
}
