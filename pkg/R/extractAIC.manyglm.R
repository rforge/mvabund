extractAIC.manyglm<- function(object, k=2, ...)
{
  n.vars = NCOL(object$y)
  ll=logLik(object)
  edf = attr(ll,"df") * n.vars
  AIC = -2*sum(ll) + k * edf
  return( c(edf,AIC) )
}
