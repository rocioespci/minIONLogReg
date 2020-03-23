\name{glm.predict.}
\alias{glm.predict.}
\title{glm.predict.}
\usage{
glm.predict.(x,slope,intercep)
}
\description{
Computes 1/(1+exp( -( slope*x + intercept ) ))
}
\arguments{
\item{x}{Numeric.}
\item{slope}{Numeric. Only one value.}
\item{intercep}{Numeric. Only one value.}
}
\details{
This function is intended to evaluate a logistic regression fit.
}
\examples{
#Compute one value
glm.predict.(x=5, slope=1, intercept=2)

#Compute several values
x <- seq(-10,10,by=0.5)
y <- glm.predict.(x=s, slope=1, intercept=2)
plot(x,y)

}

