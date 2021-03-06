% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/am_sims.R
\name{simulate_gen_assoc_recalc_PGS}
\alias{simulate_gen_assoc_recalc_PGS}
\title{Provide a quick report on the simulation of parental pairs and a child.}
\usage{
simulate_gen_assoc_recalc_PGS(
  am_corr_coeff,
  ntrait_loci,
  nloci,
  envir_var,
  pop_size,
  generations1,
  generations2,
  ntraining,
  af_r_func = runif
)
}
\arguments{
\item{am_corr_coeff}{The associative mating correlation coefficent}

\item{ntrait_loci}{The number of loci for each trait}

\item{pop_size}{The Size of the populations to simulate.  Equal numbers of males and
females of this size.}

\item{generations1}{The number of generations to run the simulation before the break.}

\item{generations2}{The number of generations after the two new polygenic scores are calculated}

\item{ntraining}{The size of the training data set}

\item{af_r_func}{A function of one variable (n) that returns a vector of allele frequencies of length n}

\item{n_loci}{The number of loci to simulate.}

\item{environ_var}{environmental variance.  A vector of length 1 or 2.}
}
\description{
This report is only for a single generation
}
