\name{meta-package}
\alias{meta-package}
\alias{meta}
\docType{package}
\title{
\packageTitle{meta}
}
\description{
Meta-Analysis for overall significance of a variable across multiple data frames
}
\details{

{Performs appropriate tests for analyzing the overall significance of a variable across multiple data frames}
\packageIndices{meta}
}
\author{
\packageAuthor{meta}

Maintainer: \packageMaintainer{meta}
}
\references{
Chang, Lun-Ching, et al. "Meta-Analysis Methods for Combining Multiple Expression Profiles: Comparisons, Statistical Characterization and an Application Guideline." BMC Bioinformatics, vol. 14, no. 1, 21 Dec. 2013, p. 368., doi:10.1186/1471-2105-14-368.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{ package }
\seealso{
}
\examples{
## generate sample data
set.seed(123)
p = 100        # total number of variables
data1 = data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(p*200), ncol=p))
data2 = data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(p*150), ncol=p))

single.variable.analysis(data1$X1, data1$group, "Mean")
oneframe.analysis(data2, "Remove")
meta.analysis(data1, data2, "Fisher")
}
