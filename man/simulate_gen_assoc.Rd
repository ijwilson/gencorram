% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/am_sims.R
\name{simulate_gen_assoc}
\alias{simulate_gen_assoc}
\title{Provide a quick report on the simulation of parental pairs and a child.}
\usage{
simulate_gen_assoc(
  am_corr_coeff,
  ntrait_loci,
  nloci,
  envir_var,
  pop_size,
  generations,
  ntraining,
  af_r_func = runif
)
}
\arguments{
\item{am_corr_coeff}{The associative mating correlation coefficent}

\item{ntrait_loci}{The number of loci for each trait}

\item{pop_size}{The Size of the populations to simulate.  Equal numbers of males and
females of this size.}

\item{generations}{The number of generations to run the simulation.}

\item{ntraining}{The size of the training data set}

\item{af_r_func}{A function of one variable (n) that returns a vector of allele frequencies of length n}

\item{n_loci}{The number of loci to simulate.}

\item{environ_var}{environmental variance.  A vector of length 1 or 2.}
}
\description{
This report is only for a single generation
}
