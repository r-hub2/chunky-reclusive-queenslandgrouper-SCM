% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.gam.R
\name{summary.scm}
\alias{summary.scm}
\title{Processing the summary to be printed}
\usage{
\method{summary}{scm}(object, intercept = FALSE, print = TRUE, ...)
}
\arguments{
\item{object}{an object of class scm}

\item{intercept}{set to TRUE if you want print the intercept}

\item{print}{print the summary}

\item{...}{further arguments to be passed to gam}
}
\value{
The summary is printed
}
\description{
This function allows to modify the summary.gam function, in order to have a bespoke summary
}
