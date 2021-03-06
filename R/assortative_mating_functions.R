#' Sample the starting genotypes
#'
#' @param num_inds  Number of individuals
#' @param allele_freq Allele frequency for each locus
#'
#' @return A matrix of allele counts, individuals in columns, loci in rows.
#' @export
#' @examples
#' a <- sample_initial_genotypes(100, rep(0.5, 20))
#' allele_freq <- runif(20, 0.01, 0.99)
#' b <- sample_initial_genotypes(100, allele_freq))
#' d <- rowMeans(b)/2
#' plot(allele_freq, d)
#'
sample_initial_genotypes <- function(num_inds, allele_freq) {
  r <- replicate(num_inds,
            rbinom(length(allele_freq), 2, allele_freq)
            )
  class(r) <- "genotypeMatrix"
  return(r)
}
#'
#' Get the genotype effects for each locus for two traits.
#'
#' The genotype effects consist of the effect for each locus, which is drawn from a Normal(0,1)
#' distribution for each locus and the  genetic standard deviation of that trait which is the
#' standard deviation for each trait normalised by the allele frequency.
#'
#' @param loci A vector with the loci affecting the trait.
#' @param allele_freq The allele frequency of sampled loci.
#' A vector at least as long as max loci specified.
#' @param env_var The environmental variance for the trait.
#' @return A list with elements of the effect and the genetic standard deviation on the trait of class effects
#' @export
generate_effects <- function(loci, allele_freq, env_var) {
  effects <- rnorm(length(loci))
  genetic_sd <- sqrt(sum((effects^2) * 2 * allele_freq[loci] * (1 - allele_freq[loci])))
  ee <- list(
      loci=loci, effects = effects, gen_sd = genetic_sd, environ_var=env_var
  )
  class(ee) <- "effects"
  return(ee)
}
#'
#' Get correlated effects for a pair of traits.
#'
#' The genotype effects consist of the effect for each locus, which is drawn from a Normal(0,1)
#' distribution for each locus and the genetic standard deviation of that trait which is the
#' standard deviation for each trait normalised by the allele frequency.
#'
#' @param loci A vector with the loci affecting the traits
#' @param allele_freq The allele frequency of sampled loci.
#' A vector at least as long as max loci specified.
#' @param corr_coeff The correlation coefficient between traits
#' @param env_var The environmental variance.  Either one value shared by both traits or
#' a vector of length 2.
#' @return A a list with effects_A and effects_B
#' @export
#' @examples 
#' f <- runif(1000, 0.05, 0.5)
#' pleio <- generate_bio_pleiotrophy_effects(1:1000, f, 0.4, env_var = c(1,1))
#' plot(pleio$effects_A$effects, pleio$effects_B$effects, xlab="Effect of locus on trait A", ylab="Effect of locus on trait B")
#' cor(pleio$effects_A$effects, pleio$effects_B$effects)
#' g <- sample_initial_genotypes(1000, f)
#' ## generate phenotypes from the model
#' phenoA <- generate_phenotype(pleio$effects_A, g)
#' phenoB <- generate_phenotype(pleio$effects_B, g)
#' ## correlation between phenotypes
#' cor(phenoA$pheno, phenoB$pheno)
#' ## correlation between additive genetic effects
#' cor(phenoA$additive_genetic, phenoB$additive_genetic)

generate_bio_pleiotrophy_effects <- function(loci, allele_freq, corr_coeff, env_var) {
  if (length(env_var) == 1)
    env_var <- c(env_var, env_var)
  
  effects <- mvtnorm::rmvnorm(length(loci), c(0, 0), 
                              matrix(c(1, corr_coeff, corr_coeff, 1), ncol=2))
  effects1 <- list(loci=loci, 
                   effects = effects[,1], 
                   gen_sd = sqrt(sum((effects[,1]^2) * 2 * allele_freq[loci] * (1 - allele_freq[loci]))),
                   environ_var=env_var[1]
  )
  effects2 <- list(loci=loci, 
                   effects = effects[, 2], 
                   gen_sd = sqrt(sum((effects[, 2]^2) * 2 * allele_freq[loci] * (1 - allele_freq[loci]))),
                   environ_var=env_var[2]
  )
  class(effects1) <- "effects"
  class(effects2) <- "effects"
  return(list(effects_A = effects1, effects_B = effects2))
}

#'
#' Get the phenotype and the additive genetic part of the trait
#'
#' Gets the phenotype for each individual for the a genetic
#' effect which comprises the additive genetic effect, and the
#' environmental effect.e
#'
#' @param effects The effect size object
#' @param genotype The genotypes with individuals in columns and sites in rows.
#' @export
#' @return A data frame with the phenotype and additive genetic variance.
#'
generate_phenotype <-  function(eff, genotype) {
  additive_genetic <- colSums(genotype[eff$loci, ] * eff$effects)

  additive_genetic <- additive_genetic / eff$gen_sd
  pheno <- additive_genetic +
    rnorm(length(additive_genetic), mean = 0.0, sd = sqrt(eff$environ_var))
  pheno <- scale(pheno, scale = FALSE) - mean(pheno)   # centre
  d <- data.frame(pheno = pheno, additive_genetic = additive_genetic)
  class(d) <- "phenotype"
  return(d)
}

#' Generate a gamete based on the parental genotype
#'
#' @param parental_genotype The genotype (as allele count) of the parent
#' @export
#' @return A vector with the haplotype (as allele count)

make_a_gamete <- function(parental_genotype) {
  u <- parental_genotype==1
  hap <- parental_genotype/2        ## works for 0 and 2
  hap[u] <- rbinom(sum(u), 1, 0.5)  ##
  return(hap)
}
#'
#' Get the next generation of genotypes
#'
#' @param mother_genotypes loci affecting trait A
#' @param father_genotypes loci affecting trait B
#' @export
#' @examples
#'
#'
next_generation <- function(mother_genotypes, father_genotypes) {
  dad_gamete <- apply(father_genotypes, 2, make_a_gamete)
  mum_gamete <- apply(mother_genotypes, 2, make_a_gamete)
  gen <- dad_gamete + mum_gamete
  class(gen) <- "genotype"
  return(gen)
}
#'
#' Get the next generation of Genotypes
#'
#' @param g a list containing male_genotype and female_genotype
#' @param effects_A the effects for trait A
#' @param effects_B the effects for trait B
#' @param corr_coeff  Correlation coefficient between traits for assortative mating
#' @export
#' @examples
#'
#'
generation <- function(g, effects_A, effects_B, corr_coeff) {
  nmales <- ncol(g$male_genotype)
  nfemales <- ncol(g$female_genotype)
  female_phenotype <- generate_phenotype(effects_A, g$female_genotype)
  male_phenotype <-   generate_phenotype(effects_B, g$male_genotype)

  # Now select the dads.  This is assuming female mate choice (I believe)
  parents_selected <- assort_mating(nmales+nfemales,
                                    female_phenotype$pheno,
                                      male_phenotype$pheno,
                                      corr_coeff)

  ## Generate - no ordering here so first half males and second half females.

  child_genotype <- next_generation(
    g$female_genotype[, parents_selected$mothers_index],
    g$male_genotype[, parents_selected$fathers_index])

  list(female_genotype=child_genotype[,1:nfemales],
       male_genotype = child_genotype[, (nfemales+1):(nmales+nfemales)])
}
#'
#' Get the next generation of Genotypes
#'
#' @param nchildren The number of children to generate
#' @param g a list containing male_genotype and female_genotype
#' @param effects_A the effects for trait A
#' @param effects_B the effects for trait B
#' @param corr_coeff  Correlation coefficient between traits for assortative mating
#' @export
#' @examples
#'
#'
sexual_reproduction <- function(nchildren, g, effects_A, effects_B, corr_coeff) {
  nmales <- ncol(g$male_genotype)
  nfemales <- ncol(g$female_genotype)
  nloci <- nrow(g$male_genotype)
  female_phenotype <- generate_phenotype(effects_A, g$female_genotype)
  male_phenotype <-   generate_phenotype(effects_B, g$male_genotype)
  
  # Now select the dads.  This is assuming female mate choice (I believe)
  parents_selected <- assort_mating(nchildren,
                                    female_phenotype$pheno,
                                    male_phenotype$pheno,
                                    corr_coeff)
  # This is not finished yet but can be speeded up considerably
  # by using Rcpp
  gg <- matrix(integer(0), ncol=nchildren, nrow=nloci)
  for (i in 1:nchildren) {
    code <- g$female_genotype[, parents_selected$mothers_index] + 3*g$male_genotype[, parents_selected$fathers_index]
    
  }
  
  ## Generate - no ordering here so first half males and second half females.
  
  child_genotype <- next_generation(
    g$female_genotype[, parents_selected$mothers_index],
    g$male_genotype[, parents_selected$fathers_index])
  
  list(female_genotype=child_genotype[,1:nfemales],
       male_genotype = child_genotype[, (nfemales+1):(nmales+nfemales)])
}
#' Generate the pairs selected by mothers
#'
#' @param pop_size The number of pairs to generate.
#' @param female_phenotype A numeric vector giving the phenotypes of all females
#' @param male_phenotype A numeric phenotype giving the phenotypes of all males
#' @param corr_coeff The correlation coefficient between phenotypes for males and females
#' @export
#' @return the indices of the returned mothers and fathers
assort_mating <-   function(popsize,
                            female_phenotype,
                            male_phenotype,
                            corr_coeff) {

  n_females <- length(female_phenotype)
  n_males <- length(male_phenotype)
  nn <- mvtnorm::rmvnorm(popsize, c(0,0), matrix(c(1, corr_coeff, corr_coeff, 1), ncol=2))
  pp <- pnorm(nn)   
  
  om <- order(male_phenotype)
  of <- order(female_phenotype)
  
  data.frame(mothers_index=of[1+floor(pp[,1]*n_females)],
                    fathers_index=om[1+floor(pp[,2]*n_males)])
}

#' Plot a phenotype
#'
#' @param phen  the phenotype
#' @export
#' @return NULL, run for side effects
#' @examples
#' allele_freq <- runif(2000, 0.01, 0.2)
#' plot(generate_phenotype(eff=generate_effects(1:1000, allele_freq, 1), sample_initial_genotypes(3000, allele_freq)))
plot.phenotype <- function(phen) {
  plot(x = phen$additive_genetic, y = phen$pheno, ylab="Phenotype", xlab="Additive Genetic" )
}
#' Find the heritability for a trait
#'
#' @param phen  The Phenotype.  An object of type phenotype
#' @export
#' @return the heritability
#' @examples
#' allele_freq <- runif(2000, 0.01, 0.2)
#' heritability(generate_phenotype(eff=generate_effects(1:1000, allele_freq, 1), sample_initial_genotypes(3000, allele_freq)))
heritability <- function(phen) {
  return(var(phen$additive_genetic)/var(phen$pheno))
}


# Code rewritten from that on Joe Pickrell blog
#
# https://github.com/joepickrell/rg-post/blob/master/Assortative_mating_sims.R
