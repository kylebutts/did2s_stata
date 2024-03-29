% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/did2s.R
\name{did2s}
\alias{did2s}
\title{Calculate two-stage difference-in-differences following Gardner (2021)}
\usage{
did2s(
  data,
  yname,
  first_stage_formula,
  treat_formula,
  treat_var,
  cluster_vars = NULL
)
}
\arguments{
\item{data}{the dataframe containing all the variables}

\item{yname}{Outcome variable}

\item{first_stage_formula}{fixed effects and other covariates you want to residualize with in first stage, use i() for fixed effects., following fixest::feols.}

\item{treat_formula}{second stage, these should be the treatment indicator(s) (e.g. treatment variable or es leads/lags), use i() for factor variables, following fixest::feols.}

\item{treat_var}{a variable that = 1 if treated, = 0 otherwise}

\item{cluster_vars}{what variable to cluster standard errors}
}
\value{
list containing fixest estimate and corrected variance-covariance matrix
}
\description{
Calculate two-stage difference-in-differences following Gardner (2021)
}
\examples{
# Load Example Dataset
data("df_hom")

# Static
static <- did2s(df_hom, yname = "dep_var", first_stage_formula = "i(state) + i(year)", treat_formula = "i(treat)", treat_var = "treat", cluster_vars = "state")
summary(static$estimate, .vcov = static$adj_cov)

# Event-Study
es <- did2s(df_hom, yname = "dep_var", first_stage_formula = "i(state) + i(year)", treat_formula = "i(rel_year)", treat_var = "treat", cluster_vars = "state")
summary(es$estimate, .vcov = es$adj_cov)

}
