\name{meta.analysis}
\alias{meta.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Meta-Analysis for overall significance of a variable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs appropriate tests to examine the overall significance of a categorical variable across different data frames
}
\usage{
meta.analysis(data1 = NULL, data2 = NULL, data3 = NULL, data4 = NULL, data5 = NULL, method = "Fisher", dataclean = "Remove")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
%%     ~~Describe \code{data1} here~~
(required) a data frame where first column consists of group numbers and all other columns are different variables. For each variable, there must be at least two unique groups. 
}
  \item{data2}{
%%     ~~Describe \code{data2} here~~
(required) a data frame with the same format as data1. Must have the same number of columns as data1.
}
  \item{data3}{
%%     ~~Describe \code{data3} here~~
(optional) a data frame with the same format as data1. If included, must have the same number of columns as the other data frames.
}
  \item{data4}{
%%     ~~Describe \code{data4} here~~
(optional) a data frame with the same format as data1. If included, must have the same number of columns as the other data frames.
}
  \item{data5}{
%%     ~~Describe \code{data5} here~~
(optional) a data frame with the same format as data1. If included, must have the same number of columns as the other data frames.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
(required) a method to be used for pooling the p-values from different data frames. Must be one of {"Fisher", "Stouffer", "Min"} or {"Max"}. Defaults to {"Fisher"}. See reference for more details.
}
  \item{dataclean}{
%%     ~~Describe \code{dataclean} here~~
a function which indicates what should happen when the input data contains {NA}s, must be one of {"Remove", "Mean"} or {"Median"}. Defaults to {"Remove"}. See details for more description.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
First column of each data frame must indicate group membership which takes at least two unique values (e.g. 1 for group1 and 2 for group2). Function will not run if the first column does not contain integer values for group memberships.

All data frames must have identical number of columns where the first column indicates group membership and columns 2 to the last column represent the variables.

Each data frame must have consistent number of rows for all columns while the number of rows may be different for different data frames. In other words, the number of samples must be equal in a single data frame but may be different across different data frames.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A summary table of the meta-analysis for the given data frames. Includes the data used, method used for pooling p-values and the test statistics and p-values for each variable with the alternative hypothesis stated at the end. You can also access the following components. See examples for more details.
  \item{statistics}{the test statistics for each column (variable) for mean comparison.}
  \item{p.value}{the p-value for the test.}
  \item{data}{the data used in the test. If there are {NA}s in the data, it will be cleaned with the method specified by the user.}
  \item{method}{the method used for pooling the p-values from different data frames. This is specified by the user.}
}
\references{
%% ~put references to the literature/web site here ~
Chang, Lun-Ching, et al. "Meta-Analysis Methods for Combining Multiple Expression Profiles: Comparisons, Statistical Characterization and an Application Guideline." BMC Bioinformatics, vol. 14, no. 1, 21 Dec. 2013, p. 368., doi:10.1186/1471-2105-14-368.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
You may print results with {print.meta()} or simply with {print}. See examples for more details.
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
print(a)
print.meta(a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
