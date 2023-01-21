\name{summary.meta}
\alias{summary.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summary of the Meta-analysis for significance of a variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Summarizes the results from the meta-analysis performed on different data frames
}
\usage{
summary.meta(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Results from the {meta.analysis()} function. See examples for more details.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A summary of the meta-analysis for the given data frames. Includes similar components as those in {meta.analysis()} function but much clearer with more useful statements at the end. You can also access the following components. See examples for more details.
  \item{data}{the data used in the test. If there are {NA}s in the data, it will be cleaned with the method specified by the user.}
  \item{method}{the method used for pooling the p-values from different data frames. This is specified by the user.}
  \item{result}{a data frame with all test statistics and p-value for each column (variable).}
  \item{samplesize}{the sample size for each data frame.}
  \item{sig1}{a data frame with all variables with p-value less than {0.01}. Also shows the test statistics and the p-values.}
  \item{sig2}{a data frame with all variables with p-value less than {0.05}. Also shows the test statistics and the p-values.}
  \item{sig3}{a data frame with all variables with p-value less than {0.10}. Also shows the test statistics and the p-values.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
You may print results with {print.summary.meta()} or simply with {print}. See examples for more details.
The printed result contains special symbols to indicate considerable p-values. {" *** "} for p-values less than {0.001}, {" ** "} for p-values between {0.001} and {0.01}, {" * "} for p-values between {0.01} and {0.05}, {" . "} for p-values between {0.05} and {0.1} and no symbol for p-values above {0.1}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## generate sample data
set.seed(123)
p = 100        # total number of variables
data1 = data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(p*200), ncol=p))
data2 = data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(p*150), ncol=p))

a = meta.analysis(data1, data2)
summary.meta(a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
