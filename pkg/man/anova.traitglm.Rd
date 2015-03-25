\name{anova.traitglm}
\alias{anova.traitglm}
\title{Analysis of deviance for a fourth corner model}
\description{
return an analysis of deviance from a fourth corner model, as computed using \code{\link{traitglm}}, to test for an environment-by-trait interaction. Slowly! This function works via \code{\link{anova.manyglm}}, which uses row-resampling for inference, and it only applies to \code{\link{traitglm}} objects that have been fitted using the (default) \code{\link{manyglm}} function.
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\usage{
\method{anova}{traitglm}(object, \dots, nBoot=100, show.time="all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted object of class \code{traitglm} and class \code{manyglm}. 
}
  \item{\dots}{
other generic \code{anova} or \code{anova.manyglm} arguments. Note however that additional \code{traitglm} objects CAN NOT be specified here, at the moment.
}
  \item{nBoot}{
The number of bootstrap iterations. Default is 100 (NOTE: increase this for later runs!)
}
\item{show.time}{Whether to display timing information for the resampling procedure: this is advisable, as resampling fourth corner models many times can take a while. The default value \code{"all"} shows all timing information, \code{"total"} shows only the overall time taken for the tests, and \code{"none"} shows none. }
}

\details{
The anova.traitglm function returns a table summarising the statistical significance of the fourth corner terms in a model, that is, the interaction between environment and traits in predicting abundance across taxa and sites. Only one fitted \code{\link{traitglm}} model can be specified, and all interaction terms from that model are simultaneously tested. The \code{\link{traitglm}} model must be fitted using the \code{\link{manyglm}} function, which is its default behaviour, in order to access the \code{\link{anova.manyglm}}, which does most of the work.

See \code{\link{anova.manyglm}} for details on how resampling is done, arguments can be passed to this function to change the test statistic (via \code{test}), the resampling method (via \code{resamp}). Because \code{\link{traitglm}} models are fitted by first vectorising the data into a univariate model, arguments such as \code{p.uni} and \code{cor.type} are redundant.

\code{\link{traitglm}} fits a single model to abundances across all sites and taxa at the same time, meaning the vector of abundances is typically pretty long, and the design matrix explaining how abundance varies across taxa and sites is typically pretty large. So resampling can take yonks. Hence the default number of resamples has been set at \code{nBoot=100}, but please consider increasing this once you have a sense for how long it will take to run (scales roughly linearly with \code{nBoot}).
}
\value{
A list of values as returned by \code{\link{anova.manyglm}}, of which the most relevant element is \code{table} (the analysis of deviance table).
}

\references{
Warton DI, Shipley B & Hastie T (in press) CATS regression - a model-based approach to studying trait-based community assembly, Methods in Ecology and Evolution
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link{traitglm}}
}

 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(antTraits)

# we'll fit a small fourth corner model, to a subset of the antTraits data.
# first select only species present in at least 25% of plots:
abSum = apply(antTraits$abund>0,2,mean)
ab = antTraits$abund[,abSum>0.25]
tr = antTraits$traits[abSum>0.25,]

# now fit the fourth corner model, only as a function of a couple of traits and env variables:
ft=traitglm(ab,antTraits$env[,1:3],data.frame(tr$Weber,tr$Femur))
anova(ft)

}
\keyword{models}
\keyword{regression}
