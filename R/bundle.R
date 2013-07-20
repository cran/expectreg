bundle <- function(B,DD,yy,pp,lambda,smooth,nb,center,constmat,types)
###### expectile regression according to eilers, schnabel
# parameters:
# formula - vector of responses ~ f(vector of independent,type="which base to use") + ...
# smooth - if smoothing with schall's algorithm, asymmetric cross validation or no smoothing shall be done
# lambda - smoothing penalty, important if no smoothing is done
{  
  nterms = length(nb)
  m = length(yy)
  np = length(pp)
  
  lala <- matrix(lambda, nrow=nterms, ncol=2, dimnames=list(1:nterms,c("mean","residual")))
  vector.a.ma.schall <- matrix(NA, nrow=sum(nb)+1*center,ncol=np)
  b=0
  cc = 0

    if(smooth == "schall")
    {
    	sch = schall(yy,B,0.5,DD,nb,lala[,1],constmat,center,types)
       lala[,1] = sch[[2]]
       vector.a.ma.schall[,1] <- mean.coefficients <- sch[[1]]
       diag.hat = sch[[3]]
       
      # dc = 1
      # it = 1
      # while(dc >= 0.01 && it < 100)# || dw != 0)
      # {
        # aa <- asyregpen.lsfit(yy, B, 0.5, lala[,1], DD, nb,constmat)
        # mean.coefficients <- aa$a 
        # diag.hat = aa$diag.hat.ma 
              
        # sig.med <- vector()
        # tau.med <- vector()
        
        # lmed <- lala[,1]
        
        # for(i in 1:nterms)
        # {
          # partbasis = (sum(nb[0:(i-1)])+1):(sum(nb[0:i]))
          # if(center)
          # {
            # partB = B[,-1,drop=FALSE][,partbasis,drop=FALSE]
            # partDD = DD[,-1,drop=FALSE][-1,,drop=FALSE][,partbasis,drop=FALSE]
            # partaa = aa$a[-1][partbasis]
          # }
          # else
          # {
            # partB = B[,partbasis,drop=FALSE]
            # partDD = DD[,partbasis,drop=FALSE]
            # partaa = aa$a[partbasis]
          # }
        
          # v <- partDD %*% partaa
          # z <- aa$fitted
       
          # H = solve(t(partB)%*%(partB) + lala[i,1]*t(partDD)%*%partDD)
          # ##H = B%*%H%*%t(B)
          # ##H = diag(H)*aa$weight
          # H = apply(partB,1,function(x){t(x)%*%H%*%x})
        
          # sig.med[i] <- sum(0.5* (yy - z) ^ 2,na.rm=TRUE) / (m - sum(aa$diag.hat.ma,na.rm=TRUE))
          # tau.med[i] <- sum(v ^ 2,na.rm=TRUE) / sum(H,na.rm=TRUE) + 1e-06
          
          # lala[i,1] <- max(sig.med[i] / tau.med[i], 1e-10,na.rm=TRUE)
        # }
        
        # dc <- max(abs(log10(lmed + 1e-6)-log10(lala[,1] + 1e-6)))
        # it = it + 1
      # }
      # if(it == 100)
        # warning("Schall algorithm did not converge. Stopping after 100 iterations.")
        
      
      #residuals = yy- aa$fitted
      residuals = yy - sch[[4]]
      
      dc = 1
      it = 1
      while(dc >= 0.01 && it < 100 && any(pp != 0.5))# || dw != 0)
      {

        b <- rep(1,ncol(B))
        cc <- pp - 0.5
        if(any(cc != 0))
        for(i in 1:20){
          mo <- fitampllsfit(residuals, B, b, pp, cc, DD, lala[,2], nb)
          b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
          c0 <- cc
          cc <- fitasy(residuals, B, b, pp, cc)
          dc <- max(abs(cc - c0))
          if (dc < 1e-6) break  
        }  

        for(q in 1:np)
        {
          vector.a.ma.schall[,q] = mean.coefficients + cc[q]*b
        }
        
        sig.res <- vector()
        tau.res <- vector()
        
        lres <- lala[,2]
        
        for(i in 1:nterms)
        {
          partbasis = (sum(nb[0:(i-1)])+1):(sum(nb[0:i]))
          if(center)
          {
            partB = B[,-1][,partbasis,drop=FALSE]
            partDD = DD[,-1][-1,][,partbasis,drop=FALSE]
            partb = b[-1][partbasis]
          }
          else
          {
            partB = B[,partbasis,drop=FALSE]
            partDD = DD[,partbasis,drop=FALSE]
            partb = b[partbasis]
          }          
          
          #partb = b[-1][partbasis]
          v = partDD %*% partb
          z = B %*% b
          
          H = solve(t(partB)%*%(partB) + lala[i,2]*t(partDD)%*%partDD)
          ##H = B%*%H%*%t(B)
          ##H = diag(H)*aa$weight
          H = apply(partB,1,function(x){t(x)%*%H%*%x})
          
          sig.res[i] <- sum(0.5*(residuals - z) ^ 2,na.rm=TRUE) / (m - sum(mo$hat.ma,na.rm=TRUE))
          #sig.res[i] = mo$sig
          tau.res[i] <- sum(v ^ 2,na.rm=TRUE) / sum(H,na.rm=TRUE) + 1e-06

          lala[i,2] <- max(sig.res[i] / tau.res[i], 1e-10,na.rm=TRUE)
        }
        
        dc <- max(abs(log10(lres + 1e-6) - log10(lala[,2] + 1e-6)))
        it = it + 1
      }
      if(it == 100)
        warning("Schall algorithm did not converge. Stopping after 100 iterations.")
    }
    else if(smooth == "gcv")
    {
      #acv.min = nlm(acv.bundle,p=lala,yy=yy,B=B,pp=pp,DD=DD,nb=nb,ndigit=8,iterlim=50,gradtol=0.0001)
      acv.min = nlminb(start=lala[,1],objective=acv,yy=yy,B=B,quantile=0.5,DD=DD,nb=nb, constmat=constmat,lower=0,upper=10000)
      # print(acv.min$estimate) 
      #min.lambda = matrix(abs(acv.min$estimate),ncol=2)
        
      #aa <- asyregpen.lsfit(yy, B, 0.5, min.lambda[,1], DD, nb)
      aa <- asyregpen.lsfit(yy, B, 0.5, abs(acv.min$par), DD, nb,constmat)
      
      mean.coefficients <- aa$a  
      lala[,1] <- abs(acv.min$par)
      diag.hat = aa$diag.hat.ma 

      residuals = yy-B%*%mean.coefficients
      
      constmat[,] = 0

      #acv.min = nlm(acv,p=lala[,2],yy=residuals,B=B,quantile=0.5,DD=DD,nb=nb,constmat=constmat,ndigit=8,iterlim=50,gradtol=0.0001)
            acv.min = nlminb(start=lala[,2],objective=acv,yy=residuals,B=B,quantile=0.5,DD=DD,nb=nb, constmat=constmat,lower=0,upper=10000)

      #print(acv.min$estimate)
      lala[,2] <- abs(acv.min$par)
      
      b <- rep(1,ncol(B))
      cc <- pp - 0.5
      if(any(cc != 0))
      for(i in 1:20){
        #mo <- fitampllsfit(residuals, B, b, pp, cc, DD, min.lambda[,2], nb)
        mo <- fitampllsfit(residuals, B, b, pp, cc, DD, abs(acv.min$par), nb)
        b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
        c0 <- cc
        cc <- fitasy(residuals, B, b, pp, cc)
        dc <- max(abs(cc - c0))
        if (dc < 1e-6) break  
      }  
  
      if(any(pp != 0.5))
      for(q in 1:np)
      {
        vector.a.ma.schall[,q] = mean.coefficients + cc[q]*b
      }
      else
        vector.a.ma.schall[,1] = mean.coefficients

    }
    else
    {

      aa <- asyregpen.lsfit(yy, B, 0.5, lala[,1], DD, nb,constmat)
      mean.coefficients <- aa$a   
      diag.hat = aa$diag.hat.ma 

      residuals = yy-B%*%mean.coefficients

      b <- rep(1,ncol(B))
      cc <- pp - 0.5
      if(any(cc != 0))
      for(i in 1:20){
        mo <- fitampllsfit(residuals, B, b, pp, cc, DD, lala[,2], nb)
        b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
        c0 <- cc
        cc <- fitasy(residuals, B, b, pp, cc)
        dc <- max(abs(cc - c0))
        if (dc < 1e-6) break  
      }  
  
      for(q in 1:np)
      {
        vector.a.ma.schall[,q] = mean.coefficients + cc[q]*b
      }
  
  }
  
  diag.hat = matrix(diag.hat,nrow=length(diag.hat),ncol=np)

  return(list(vector.a.ma.schall,lala,diag.hat,mean.coefficients,b,cc))
}




fitampllsfit <- function(y, B, b, p, cc, DD, lambda, nb)
{
  # Fit amplitude
  
  # Initialize
  m <- nrow(B)
  n <- ncol(B)
  np <- length(p)
  a <- B %*% b
  BtB <- 0
  Bty <- 0
  
  w = matrix(NA,nrow=length(y),ncol=np)
  
  # Add inner products for all p
  for (j in 1:np){
    #w[,j] <- ifelse(y >= cc[j]* a, p[j], 1 - p[j])
    w[,j] = p[j]
    w[!(y >= cc[j]*a),j] = 1 - p[j]

    WB <- cc[j] * as.vector(w[,j]) * B
    BtB <- BtB + cc[j] * t(WB) %*% B
    Bty <- Bty + t(WB) %*% y
  }

  lambda = c(rep(0,times=n - sum(nb)),rep(lambda,times=nb))

  # Solve penalized equations
  P <- sqrt(lambda) * t(DD) %*% DD

  model <- lsfit(x=BtB+P,y=Bty, intercept=FALSE)

  sigma = .5*sum((Bty - (BtB+P) %*% model$coef)^2,na.rm=TRUE) / (length(Bty) - sum(hat(model$qr)[1:length(Bty)]))

  return(list(b=model$coef, hatma = hat(model$qr)[1:length(Bty)], weight=w,sig=sigma))
}




fitasy <- function(y, B, b, p, c0)
{
  a <- B %*% b
  ccr <- 0 * p
  for(j in 1:length(p)){
    w <- ifelse(y >= c0[j]*a, p[j], 1 - p[j])
    ccr[j] <- sum(w * a * y) / sum(w * a * a)
  }
  return(ccr)
}


acv.bundle <- function(penalty,yy,B,pp,DD,nb)
# asymmetric cross validation
# computes the acv score for the smoothing of the regression
# score has to be minimized dependant on parameter "penalty"
# therefore a grid search can be applied to this function
# parameters:
# penalty - smoothing parameter lambda
# yy - vector of responses
# B - basis for the approximation
# pp - quantile
# DD - penalization matrix
# nb - vector with number of basis elements for each term
{
  penalty = matrix(abs(penalty),ncol=2)
  aa <- asyregpen.lsfit(yy, B, 0.5, penalty[,1], DD, nb)

  residuals = yy-B%*%aa$a

  b <- rep(1,ncol(B))
  cc <- pp - 0.5
  for(i in 1:20){
    mo <- fitampllsfit(residuals, B, b, pp, cc, DD, penalty[,2], nb)
    b <- mo$b/((sum(mo$b^2)/length(mo$b))^(1/2))
    c0 <- cc
    cc <- fitasy(residuals, B, b, pp, cc)
    dc <- max(abs(cc - c0))
    if (dc < 1e-6) break  
  }  
  
  coef.vector = matrix(NA, nrow=sum(nb)+1,ncol=length(pp))
  for(q in 1:length(pp))
  {
    coef.vector[,q] = aa$a + cc[q]*b
  }

  #aa <- asyregpen.lsfit(yy, B, quantile, abs(penalty), DD, nb)

  #H = solve(t(B)%*%(aa$weight*B) + penalty*t(DD)%*%DD)
  ##H = B%*%H%*%t(B)
  ##H = diag(H)*aa$weight
  #H = apply(sqrt(aa$weight)*B,1,function(x){t(x)%*%H%*%x})

  #H = diag(diag(aa$diag.hat.ma))

  #score = 0
  #for(i in 1:length(pp))
  #{
  #  score = score + mo$weight[,i]*(yy - B%*%coef.vector[,i])^2/(1-aa$diag.hat.ma)^2
  #}
  
  score = (yy - B%*%aa$a)^2/(1-aa$diag.hat.ma)^2
  
  resid = asyregpen.lsfit(residuals, B, 0.5, penalty[,2], DD, nb)
  
  score = score + (residuals - B%*%resid$a)^2/(1-resid$diag.hat.ma)^2
    
  mean(score[which(is.finite(score))],na.rm=TRUE)
}
