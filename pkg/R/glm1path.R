glm1path = function(y, X, family="negative.binomial", lambdas=NULL, penalty = c(0, rep(1, dim(X)[2]-1)), df.max = sum(y>0),
                    n.lambda=25, lam.max=NULL, lam.min=NULL, k=log(length(y)), b.init=NA, phi.init=NA, phi.iter=1, ...)

################################################################################################
{
  
  v.inv = 1
  mu = mean(y)

  # first turn family into a proper family function, if needed

  allargs <- match.call(expand.dots = FALSE)
  dots <- allargs$...
  if( "tol" %in% names(dots) )
    tol = dots$tol
  else
    tol = c(1.e-8, .Machine$double.eps)
  if( "n.iter" %in% names(dots) )
    n.iter = dots$n.iter
  else
    n.iter = 100
  
  if(is.character(family)) 
  {
    if (family == "negbinomial" || family=="negative.binomial")
    {
      family = negative.binomial(1/tol[2])
      phi.init=NA
    }
    else
    {
      fam.fn = get(family, mode = "function", envir = parent.frame())        
      family = fam.fn()
    }
  }

  if(pmatch("Negative Binomial",family$family,nomatch=0)==1)
  {
    if(is.na(phi.init))
    {
      phi.init = family$var(1)-1 #extract overdispersion parameter by evaluating variance at mu=1 and solving for phi.
      if(phi.init>2*tol[2])
        phi.iter=n.iter*2 #if phi was specified in family argument, make sure it stays fixed in estimation.
    }
    v.inv = mu/(mu+phi.init*mu^2)
  }
  
  res = v.inv * ( y - mean(y) )
  score  = t(X)  %*% res
  max.score = max(abs(score[penalty>0]))

  # get a sensible initial estimate of phi from data, if none has been provided
  if(pmatch("Negative Binomial",family$family,nomatch=0)==1 & phi.init<=2*tol[2] )
  {
    init     = glm1(y, X, penalty*max.score, family=family, b.init=b.init, phi.init=NA, phi.iter=1, ...)
    phi.init = init$phi
    v.inv    = mu / (mu+phi.init*mu^2)
    # re-estimate max.score now phi has been updated:
    res = v.inv * ( y - mean(y) )
    score  = t(X)  %*% res
    max.score = max(abs(score[penalty>0]))
  }

  #Determine the range of lambda values to assess:
  if(is.null(lambdas))
  {
    if(is.null(lam.max))
      lam.max = max.score
    if(is.null(lam.min))
      lam.min = lam.max / 100000
    lambdas = exp( seq( log(lam.max), log(lam.min), length=n.lambda ) )
    if(lam.max!=max.score) #to add intercept model to path if not yet included
    {
      lambdas = c(max.score,lambdas)
    }
  }
  else
    lambdas=sort(lambdas,decreasing=TRUE) #make sure they are sorted in decreasing order

  n.lambda = length(lambdas)
  

################## Fit model to full dataset ########################
  ll         = rep( NA, length=n.lambda)
  phi        = ll
  df         = ll
  counter    = ll
  check      = ll
  beta       = matrix( NA, dim(X)[2],   n.lambda, dimnames = list( dimnames(X)[[2]],   lambdas ) )

  b.old    = b.init
  phi.old  = phi.init

  #estimating phi and beta jointly:
  for ( i.lambda in 1:n.lambda)
  {
    penalty.i = lambdas[i.lambda] * penalty
    out       = glm1(y, X, penalty.i, family=family, b.init=b.old, phi.init=phi.old, phi.iter=phi.iter, ...)
    b.old     = out$beta
    phi.old   = out$phi
  
    ll[i.lambda]      = out$likes[length(out$likes)]
    df[i.lambda]      = sum(abs(out$beta)>tol[1])
    beta[,i.lambda]   = out$beta
    phi[i.lambda]     = out$phi
    counter[i.lambda] = out$counter
    check[i.lambda]   = out$check
    if(df[i.lambda]>df.max)
      break    
  }
  if(i.lambda<n.lambda) #delete all the NA's if broke out of loop early:
  {
    ll=ll[1:i.lambda]
    df=df[1:i.lambda]
    beta=beta[,1:i.lambda]
    phi = phi[1:i.lambda]
    counter = counter[1:i.lambda]
    check = check[1:i.lambda]
    lambdas = lambdas[1:i.lambda]
    n.lambda=i.lambda
  }
  bics  = -2*ll + k * df

  id.use = which(bics==min(bics))[1]

  #refit best model to get all its bells and whistles
  penalty.i = lambdas[id.use] * penalty
  best = glm1(y, X, penalty.i, family=family, b.init=beta[,id.use], phi.init=phi[id.use], phi.iter=phi.iter)

  lasso.final = list(beta.best = beta[,id.use], lambda.best = lambdas[id.use], glm1.best=best, ll=ll, df=df, bics=bics, counter=counter, check=check, phi=phi, lambdas=lambdas, all.beta=beta)

}
# end function
