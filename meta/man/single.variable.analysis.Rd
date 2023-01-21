\name{single.variable.analysis}
\alias{single.variable.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Single-variable Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compares group mean values for a given single variable
}
\usage{
single.variable.analysis(x, group, dataclean = "Remove")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a list of numerical measurements for a single variable. Measurements must come from at least two different groups.
}
  \item{group}{
%%     ~~Describe \code{group} here~~
a list of group membership (i.e. 1 for group1 and 2 for group2), must have at least two unique groups and group names must be integers.
}
  \item{dataclean}{
%%     ~~Describe \code{dataclean} here~~
a function which indicates what should happen when the input data contains {NA}s, must be one of {"Remove", "Mean"} or {"Median"}. Defaults to {"Remove"}. See details for more description.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If either {x} or {group} contains {NA}s, you can choose to remove the missing value ({"Remove"}), replace it with the mean value ({"Mean"}) or the median value ({"Median"}) of the group of the missing value.

For data with two unique groups, if both groups have normal data, {t.test()} is used. If at least one group has non-normal data, {wilcox.test()} is used.

For data with more than two groups, if all groups have normal data, {anova(lm())} is used. If at least one group has non-normal data, {kruskal.test()} is used.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A summary table of the test performed. Includes the method name, test statistic and the p-value. You can also access the following components. See examples for more details.
A list containing the following components:
  \item{method}{the method used for mean comparison.}
  \item{statistic}{the test statistic of the chosen method.}
  \item{p.value}{the p-value for the test.}
  \item{data}{the data used in the test. If there are {NA}s in the data, it will be cleaned with the method specified by the user.}
  \item{model}{the final summary table of the test performed.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
You may print results with {print.singlevariable()} or simply with {print}. See examples for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{t.test}}
  \code{\link{wilcox.test}}
  \code{\link{anova}}
  \code{\link{lm}}
  \code{\link{kruskal.test}}
}
\examples{
## generate sample data
set.seed(123)
p = 100        # total number of variables
data1 = data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(p*200), ncol=p))

a = single.variable.analysis(data1$X1, data1$group)
print(a)
print.singlevariable(a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
