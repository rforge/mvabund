traitglm = function( L, R, Q=NULL, family="negative.binomial", method="manyglm", ...  )
{

#subfunctions get.design and get.polys defined below.
  
  # extract any arguments that work with cv.glm1path and save separately so they are not passed to glm1path.
  allargs <- match.call(expand.dots = FALSE)
  dots <- allargs$...
  if( "best" %in% names(dots) )
    best <- dots$best    
  else
    best = "1se"
  if( "plot" %in% names(dots) )
    plot <- dots$plot    
  else
    plot=TRUE
  if( "prop.test" %in% names(dots) )
    prop.test <- dots$prop.test
  else
    prop.test = 0.2
  if( "n.split" %in% names(dots) )
    n.split <- dots$n.split
  else
    n.split=10
  if( "seed" %in% names(dots) )
    seed <- dots$seed
  else
    seed=NULL
  if( "show.progress" %in% names(dots) )
    show.progress <- dots$show.progress
  else
    show.progress = FALSE
  if( "get.fourth" %in% names(dots) )
    get.fourth <- dots$get.fourth
  else
    get.fourth = TRUE
  
  deactive <- c("best", "plot", "prop.test", "n.split", "seed", "show.progress", "get.fourth") 
  deactivate <- (1:length(dots))[names(dots) %in% deactive ]  
  for (i in length(deactivate):1) 
    dots[ deactivate[i] ]<-NULL
  
  dots <- lapply( dots, eval, parent.frame() )
  
  
    n.sites = dim(L)[1] #use R.des gives number of sites for prediction / model fitting
    n.spp   = dim(L)[2]

    # get standardised R, Q, orthogonal polys, and poly coeffs
    R.des = get.polys(R)
    
    if(is.null(Q))
    {
        cat(paste("No traits matrix entered, so will fit SDMs with different env response for each spp","\n"))
        Q.des = list(X=data.frame(names(L)),X.squ=NULL,var.type="factor")
    }
    else
        Q.des = get.polys(Q)

    any.penalty = method=="cv.glm1path" || method=="glm1path"
    # get mega-matrix of design for regression against vectorised l
    X.des = get.design( R.des, Q.des, names(L), spp.penalty=TRUE, any.penalty=any.penalty, get.fourth=get.fourth )
    X = X.des$X
    Lm <- as.matrix(L)
    l <- as.vector(Lm)
    rm(Lm)

#    require(Matrix)
#    require(MASS)

    if( method=="cv.glm1path" || method=="glm1path" )
    {
      block = factor(rep(1:n.sites,n.spp))
      ft = do.call(glm1path,c(list(y=l, X=X, family=family, penalty = X.des$penalty), k=log(n.sites), dots))
      if( method=="cv.glm1path" )
        ft = do.call(cv.glm1path,c(list(object=ft, block=block, best=best, plot=plot, prop.test=prop.test, n.split=n.split, 
                                      seed=seed, show.progress=show.progress), dots))
      id.use = which(ft$lambdas==ft$lambda)
      ft$deviance = -2*ft$logL[id.use]
      ft$phi = ft$glm1$phi
      if(ft$df[1]==1)
        null.deviance = ft$logL[1]
#      ft$family=ft$glm1$family
    }
    else
      ft = do.call( method, c(list(formula=l~., family=family, data=data.frame(X)), dots) )

    # report fourth corner matrix for "best" model
    ft$fourth.corner = matrix( coef(ft)[X.des$is.4th.corner], length(X.des$names.Q), length(X.des$names.R) )
    dimnames(ft$fourth.corner)[[1]] = X.des$names.Q
    dimnames(ft$fourth.corner)[[2]] = X.des$names.R
#    ft$fourth.corner = Matrix(fourth.corner, sparse=T)
    # cat("\n")
    # cat("Fourth corner matrix:    ")
    # cat("\n")
    # print(round(ft$fourth.corner,dec.pl))


    ### Plot results

    # do a lattice "levelplot" of the coefficients in the matrix of fourth corner interactions
#    if( fourthPlot == TRUE & get.fourth == TRUE )
#    {
#        library(lattice)
#        a        = max( abs(ft$fourth.corner) )
#        if(a < 1.e-8)
#            warning("No fourth corner interactions were included in final model")
#        else
#        {
#            colort   = colorRampPalette(c("blue","white","red")) 
#            plot.4th = levelplot(t(as.matrix(ft$fourth.corner)), xlab="Environmental Variables", ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100), scales = list( x= list(rot = 45)))
#            print(plot.4th)
#        }
#    }
    
    ft$R.des = R.des
    ft$Q.des = Q.des
    ft$any.penalty = any.penalty
    ft$L = L
    ft$scaling = X.des$scaling
    ft$call=match.call()

    class(ft)=c("traitglm",class(ft))
    return( ft )
    
}



################ get.polys for getting orthogonal polynomials ###################
     
get.polys = function( X, X.des.train=NULL)
{
# get.polys will take a matrix of env vars (or trait vars), and standardise the quant ones
# as well as return orthogonal poly's. Importantly, if training matrices are given as input,
# these will be used in matrix construction.

    n.sites = dim(X)[1]
    n.params  = dim(X)[2]
    if(is.null(X.des.train))
        n.train.sites = n.sites
    else
        n.train.sites <- dim(X.des.train$X)[1]
    if(is.null(X.des.train$var.type))
        var.type = rep("quantitative",n.params)
    else
        var.type = X.des.train$var.type
    for (i in 1:n.params)
    {
        if(is.factor(X[,i]))
        {
            n.levels    = length(levels(X[,i]))
            if(n.levels==2)
            {
                var.type[i]="binary" #treat as quantitative but don't find its square
                #change variable name to indicate what it is doing
                dimnames(X)[[2]][i]=paste(names(X)[i],levels(X[,i])[2],sep="")
                #change entry in X to numerical
                X[,i]  = as.numeric(as.factor(X[,i]))*2-3
            }
            else
            {
                var.type[i]="factor"
                contrasts(X[,i])[contrasts(X[,i])==0] = -1
            }
        }
    }

    # to return standardised values of quant vars, with coeff, where needed:
    is.quant = which(var.type=="quantitative")
    n.quant = length(is.quant)
    if( n.quant>0 )
    {
        X.squ = X[,0]
        degs = c()
        names.X.squ = c()
        poly.coefs = as.list( rep( NA, n.quant ) )
        for(i.quant in is.quant)
        {
            poly.i = poly( X[,i.quant], degree=2, coefs=X.des.train$coefs[[i.quant]] )
            X.squ = cbind( X.squ, poly.i )
            degs  = c( degs, attr(poly.i, "degree") )
            poly.coefs[[ i.quant ]] = attr(poly.i, "coefs")   
            names.X.squ = c( names.X.squ, dimnames(X)[[2]][i.quant], paste( dimnames(X)[[2]][i.quant], ".squ", sep="") )
        }
        X.squ = X.squ * sqrt(n.train.sites)
        dimnames(X.squ)[[2]] = names.X.squ
        X[,var.type=="quantitative"] = X.squ[,degs==1]
        #get rid of the linear terms:
        X.squ = X.squ[,degs==2]
    }
    else
    {
        X.squ=NULL
        poly.coefs=NULL
    }
    # to return orthogonal poly values of quant vars (with coeff):
    return( list( X=X, X.squ=X.squ, var.type=var.type, coefs=poly.coefs ) )
}



################ get.design for getting the design matrix ###################
get.design = function( R.des, Q.des, L.names, spp.penalty=FALSE, any.penalty=TRUE, scaling=NULL, get.fourth=TRUE )
{

# get.design will take matrices of linear env and trait terms, and orthogonal quadratic terms,
# and return a mega-matrix that can be regressed against vectorised abundance.
# also returns logical for fourth corner terms, and their row and column names

    #How many site are there? This will be used later
    n.sites = dim(R.des$X)[1] #use R.des gives number of sites for prediction / model fitting
    n.spp   = length(L.names)

    is.scaling.given = is.null(scaling)==F #remember if scaling was specified - if not, build it up
    #get spp indicator
    if(is.null(spp.penalty)==FALSE)
    {
        spp     = rep(L.names,each=n.sites)
        spp     = as.factor(spp)
        mod     = as.formula("~spp-1")
        X.spp   = model.matrix(mod)
        if(spp.penalty==FALSE) #if species not penalised need to remove species 1 indicator
          X.spp   = X.spp[,-1]
        X.spp[X.spp==0] = -1
        if(is.scaling.given==F) #rescale X.spp so ALL variables have variance 1.
        {
          scaling = list()
          X.spp = scale(X.spp)
          scaling$spp$center = attr(X.spp,"scaled:center")
          scaling$spp$scale  = attr(X.spp,"scaled:scale")      
        }
        if(is.scaling.given)
          X.spp = scale(X.spp,center=scaling$spp$center, scale=scaling$spp$scale)
        X.spp   = cbind(1,X.spp) #add intercept since my code doesn't do that automatically
    }
    else
    {
        X.spp = as.matrix( rep(1,n.sites*n.spp) )
        if(is.scaling.given==F) #if required, initiate scaling object (to be built later) 
          scaling = list() 
    }

    # R terms
    X.R     = X.spp[,0]
    if( is.null(R.des$X.squ) )
    {
        R.small = R.des$X
        var.type= R.des$var.type
        is.lin.small = rep( TRUE,NCOL(R.des$X) )
    }
    else
    {
        R.small = cbind( R.des$X, R.des$X.squ )
        var.type= c( R.des$var.type, rep("quantitative",dim(R.des$X.squ)[2]) )
        is.lin.small = c( rep( TRUE,NCOL(R.des$X) ), rep( FALSE, NCOL(R.des$X.squ) ) )
    }    
    names.R = c()
    is.lin.R= c()
    for( iR in 1:NCOL(R.small) )
    {
        R.i   = rep( R.small[,iR], times=n.spp )
        mod     = as.formula("~0+R.i")
        mm      = model.matrix(mod)
        names.i = dimnames(R.small)[[2]][iR]
        if(var.type[iR]=="factor")
        {
            names.i = paste(dimnames(R.small)[[2]][iR], levels(R.small[,iR]),sep="")
        }
        names.R = c( names.R, names.i )
        is.lin.R= c( is.lin.R, rep( is.lin.small[iR] , dim(mm)[2] ) )
        X.R     = cbind( X.R, mm )
    }
    dimnames(X.R)[[2]]=names.R
    if(is.scaling.given==F) #rescale X.R so ALL variables have variance 1.
    {
      X.R = scale(X.R)
      scaling$R$center = attr(X.R,"scaled:center")
      scaling$R$scale  = attr(X.R,"scaled:scale")      
    }
    if(is.scaling.given)
      X.R = scale(X.R,center=scaling$R$center, scale=scaling$R$scale)
    
    # Q terms
    X.Q     = X.spp[,0]
    if( is.null(Q.des$X.squ) )
    {
        Q.small = Q.des$X
        var.type= Q.des$var.type
        is.lin.small = rep( TRUE,NCOL(Q.des$X) )
    }
    else
    {
        Q.small = cbind( Q.des$X, Q.des$X.squ )
        var.type= c( Q.des$var.type, rep("quantitative",dim(Q.des$X.squ)[2]) )
        is.lin.small = c( rep( TRUE,NCOL(Q.des$X) ), rep( FALSE, NCOL(Q.des$X.squ) ) )
    }    
    names.Q = c()
    is.lin.Q= c()
    for( iQ in 1:NCOL(Q.small) )
    {
        Q.i  = rep( Q.small[,iQ], each=n.sites )
        mod     = as.formula("~0+Q.i")
        mm      = model.matrix(mod)
        names.i = dimnames(Q.small)[[2]][iQ]
        if(var.type[iQ]=="factor")
        {
            names.i = paste(dimnames(Q.small)[[2]][iQ], levels(Q.small[,iQ]),sep="")
        }
        names.Q = c( names.Q, names.i )
        is.lin.Q= c( is.lin.Q, rep( is.lin.small[iQ] , dim(mm)[2] ) )
        X.Q     = cbind( X.Q, mm )
    }
    dimnames(X.Q)[[2]]=names.Q
    if(is.scaling.given==F)
    {
      X.Q = scale(X.Q)
      scaling$Q$center = attr(X.Q,"scaled:center")
      scaling$Q$scale  = attr(X.Q,"scaled:scale")      
    }
    if(is.scaling.given)
      X.Q = scale(X.Q, center=scaling$Q$center, scale=scaling$Q$scale)
    
    if(get.fourth==TRUE)
    {
      # R*Q interaction
      X.RQ    = X.spp[,0]
      n.R = sum(is.lin.R)
      n.Q = sum(is.lin.Q)
      ref.R = rep(1:n.R,each=n.Q)
      ref.Q = rep(1:n.Q,n.R)
      X.RQ = X.R[,ref.R] * X.Q[,ref.Q]
      dimnames(X.RQ)[[2]] = paste(dimnames(X.R)[[2]][ref.R], dimnames(X.Q)[[2]][ref.Q], sep=":")
      if(is.scaling.given==F)
      {
        X.RQ = scale(X.RQ)
        scaling$RQ$center = attr(X.RQ,"scaled:center")
        scaling$RQ$scale  = attr(X.RQ,"scaled:scale")      
      }
      if(is.scaling.given)
        X.RQ = scale(X.RQ,center=scaling$RQ$center, scale=scaling$RQ$scale)
    }
    else
      X.RQ=X.R[,0]
    
    if(any.penalty)
      X             = cbind(X.spp,X.R,X.Q,X.RQ)
    else
      X             = cbind(X.spp,X.R,X.RQ) #no trait main effects if there is no penalty on trait params, covered by spp.
    
    n.X           = dim(X)[2]
    if(get.fourth==TRUE)
      is.4th.corner = c( rep(F,n.X-n.R*n.Q), rep(T,n.R*n.Q) )
    else
      is.4th.corner = rep(F,n.X)
    names.R = dimnames(X.R)[[2]]
    names.Q = dimnames(X.Q)[[2]]

    penalty = c( 0, rep(1,dim(X)[2]-1) )
    if (is.null(spp.penalty)==FALSE) #change if spp intercept terms are to be left unpenalised:
    {
      if (spp.penalty==FALSE)
        penalty  = c( rep( 0,dim(X.spp)[2] ), rep( 1, dim(X)[2]-dim(X.spp)[2] ) )
    }
    return(list(X=X, is.4th.corner=is.4th.corner, names.R=names.R[is.lin.R], names.Q=names.Q[is.lin.Q], penalty=penalty, any.penalty=any.penalty, scaling=scaling) )

}
