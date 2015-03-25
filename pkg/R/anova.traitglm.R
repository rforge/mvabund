anova.traitglm <- function(object, ..., nBoot=100, show.time="all")
{
    if("manyglm" %in% class(object) == FALSE)
      stop("Sorry, the anova function only works for traitglm objects fitted using manyglm.")

    n.sites    = dim(object$L)[1]
    n.spp      = dim(object$L)[2]
    block      = factor(rep(1:n.sites,n.spp))
    object$call$get.fourth = FALSE
    env.plus.trait = eval(object$call)
    env.times.trait=object
    an = anova.manyglm(env.plus.trait, env.times.trait, nBoot=nBoot, show.time=show.time, ...)
    
    #get rid of multivariate terms that don't apply, since it is a univariate fit:
    an$p.uni="none"
    an$cor.type="I"
    an$uni.p[is.numeric(an$uni.p)] <- NA
    row.names(an$table)=c("Main effects only", "env:trait (fourth corner)")
    return(an)
}
