\name{antTraits}
\docType{data}
\alias{antTraits}
\title{Ant data, with species traits}
\description{
Abundances of 41 epigaeic ant species across 30 sites in south-eastern Australia, with species trait and environmental data 
}
\usage{
data(antTraits)
}
\format{
A list containing three elements:
\describe{
\item{abund}{A data frame with observations at 30 different locations
 of abundances on 41 epigaeic ant species.
}
\item{env}{A data frame containing 7 environmental variables from transects at each of the 30 sites. Most varibles are percent cover (either in two 20x20m transects, or in ten 1x1 metre quadrats), but \code{feral.mammal.dung} reports the proporion out of 10 quadrats in which dung was found. Details in Gibb et al (2015).}
\item{traits}{A data frame containing 5 species traits measured for each of the 41 species. Weber's length was log-transformed, Femur length was log-transformed then regressed against log(Weber's length), to remove the effect of size.}
}
}
\sources{Data provided by Heloise Gibb.}
\references{
Gibb H, Stoklosa J, Warton, DI, Brown, AM, Andrew, NR and Cunningham, SA (2015) Does morphology predict trophic position and habitat use of ant species and assemblages? Oecologia 177, 519-531.
}
\examples{
data(antTraits)
ft = traitglm(ab,env2,tr) #to do a fourth corner analysis
}
\keyword{datasets}
