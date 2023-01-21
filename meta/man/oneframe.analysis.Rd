\name{oneframe.analysis}
\alias{oneframe.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
One Data frame Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compares group mean values for all variables in single data frame
}
\usage{
oneframe.analysis(x, dataclean = "Remove")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a data frame where first column consists of group numbers and all other columns are different variables. For each variable, there must be at least two unique groups.
}
  \item{dataclean}{
%%     ~~Describe \code{dataclean} here~~
a function which indicates what should happen when the input data contains {NA}s, must be one of {"Remove", "Mean"} or {"Median"}. Defaults to {"Remove"}. See details for more description.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is an expansion of the function {single.variable.analysis()} in this package. The {single.variable.analysis()} function is applied to all columns (except the first group membership column) to compare group mean values and returns p-values for each variable.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A summary table with method, test statistic and p-value for all variables. You can also access the following components. See examples for more details.
  \item{method}{the method used in each column (variable) for mean comparison.}
  \item{statistic}{the test statistic of the chosen method.}
  \item{p.value}{the p-value for the test.}
  \item{data}{the data used in the test. If there are {NA}s in the data, it will be cleaned with the method specified by the user.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
You may print results with {print.oneframe()} or simply with {print}. See examples for more details.
The printed result contains special symbols to indicate considerable p-values. {" *** "} for p-values less than {0.001}, {" ** "} for p-values between {0.001} and {0.01}, {" * "} for p-values between {0.01} and {0.05}, {" . "} for p-values between {0.05} and {0.1} and no symbol for p-values above {0.1}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{single.variable.analysis}}
}
\examples{
## generate sample data
set.seed(123)
p = 100        # total number of variables
data1 = data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(p*200), ncol=p))

a = oneframe.analysis(data1)
print(a)
print.oneframe(a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
