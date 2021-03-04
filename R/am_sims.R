#'
#' Generate father mother pairs
#'
#' This is really just a test that the simulations and the
#' correlation coefficients work.
#'
#' @param population_size The Size of the populations to simulate.  Equal numbers of males and
#' females of this size.
#' @param loci_A loci affecting trait A
#' @param loci_B loci affecting trait B
#' @param environ_var environmental variance
#' @param corr_coeff Correlation coefficient between phenotype in males and females
#' @param n_loci The number of loci to simulate.  Not always needed.
#' @param allele_freq The allele frequency of sampled loci.
#' Either a single value or a vector at least as long as max loci specified.
#' @export
#' @examples
#' g1 <- generate_QT_assortative_mating(population_size = 100000, loci_A = 1:100, loci_B = 101:200,
#'                                      environ_var = 0.3, corr_coeff = 0.3)
#' cor(g1$mother_phenotype$pheno, g1$father_phenotype$pheno)#
#' plot(g1$mother_phenotype$pheno, g1$father_phenotype$pheno)
#' # Simulate for loci with rarer variants
#' g2 <- generate_QT_assortative_mating(population_size = 100000, loci_A = 1:100, loci_B = 101:200,
#'                          environ_var = 0.3, corr_coeff = 0.3, allele_freq = runif(200, 0.01, 0.05))
#' cor(g2$mother_phenotype$pheno, g2$father_phenotype$pheno)
#' plot(g2$mother_phenotype$pheno, g2$father_phenotype$pheno)
#'
generate_QT_assortative_mating <- function(
  population_size,
  loci_A, loci_B,
  environ_var,  # environmental variance
  corr_coeff = 0,
  n_loci,
  allele_freq = 0.5) {

  g <- init(population_size, loci_A, loci_B, environ_var, n_loci, allele_freq)

  ## Get the female and male phenotypes in the initial generation
  female_phenotype <- generate_phenotype(g$effects_A, g$female_genotype)
  male_phenotype <-   generate_phenotype(g$effects_B, g$male_genotype)

  # Now select the mums and dads.  This is assuming female mate choice (I believe)
  parents <- assort_mating(2*population_size, female_phenotype$pheno,
                           male_phenotype$pheno, corr_coeff)

  mother_genotype  <- g$female_genotype[, parents$mothers_index]
  father_genotype <- g$male_genotype[, parents$fathers_index]

  print(head(female_phenotype))

  return(
    list(
      mother_geno = mother_genotype,
      father_geno = father_genotype,
      effects_A = g$effects_A,
      effects_B = g$effects_B,
      mother_phenotype = female_phenotype[parents$mothers_index,],
      father_phenotype = male_phenotype[parents$fathers_index,],
      allele_freq=allele_freq)
  )
}

#' Provide a quick report on the simulation of parental pairs and a child.
#'
#' This report is only for a single generation
#'
#' @param population_size The Size of the populations to simulate.  Equal numbers of males and
#' females of this size.
#' @param loci_A loci affecting trait A
#' @param loci_B loci affecting trait B
#' @param environ_var environmental variance
#' @param corr_coeff Correlation coefficient between phenotype in males and females
#' @param n_loci The number of loci to simulate.  Not always needed.
#' @param allele_freq The allele frequency of sampled loci.
#' Either a single value or a vector at least as long as max loci specified.
#' @export

QT_assortative_mating_report <- function(population_size,
                                  loci_A,       # loci affecting trait A
                                  loci_B,       # loci affecting trait B
                                  environ_var,  # environmental variance
                                  corr_coeff = 0,
                                  n_loci,
                                  allele_freq = 0.5,
                                  plots = FALSE,
                                  write_plots = FALSE) {

  g <- generate_QT_assortative_mating(population_size,
                                      loci_A, loci_B,
                                      environ_var, corr_coeff, n_loci, allele_freq)

  print(names(g))

  L <- length(allele_freq)
  cor_pheno <- cor(g$mother_phenotype$pheno, g$father_phenotype$pheno)
  cat("correlation coefficient between parental phenotypes", cor_pheno, "\n")
  cor_add_gen <- cor(g$mother_phenotype$additive_genetic, g$father_phenotype$additive_genetic)
  cat("correlation coefficient between parental additive phenotype", cor_add_gen, "\n")

  layout(t(1:2))

  plot(
    g$mother_phenotype$pheno, g$father_phenotype$pheno,
    xlab = "Female Phenotype",
    ylab = "Male Phenotype",
    col = adjustcolor("black", 0.2),
    main = paste("Corr. parental phenotype=",  corr_coeff)
  )
  abline(lm(g$father_phenotype$pheno ~ g$mother_phenotype$pheno, col = "red"))
  text(x = -2,
       y = 2,
       format(cor(g$father_phenotype$pheno, g$mother_phenotype$pheno), dig = 2),
       col = "red")

  plot(
    g$mother_phenotype$additive_genetic, g$father_phenotype$additive_genetic,
    xlab = "Female Additive Genetic Score A",
    ylab = "Male Additive Genetic Score B",
    col = adjustcolor("black", 0.2),
    main = paste("Corr. parental additive genetic score=",  corr_coeff)
  )
  abline(lm(g$father_phenotype$additive_genetic ~ g$mother_phenotype$additive_genetic, col = "red"))
  text(x = -2,
       y = 2,
       format(cor(g$father_phenotype$additive_genetic, g$mother_phenotype$additive_genetic), dig = 2),
       col = "red")



  plot(
    colSums(g$mother_geno) / L,
    colSums(g$father_geno) / L,
    xlab = "Female Mean Genotype",
    ylab = "Male Mean Genotype",
    col = adjustcolor("black", 0.2),
    main = paste(
      "Cor. parental pheno.=",
      corr_coeff,
      "Heritability",
      format(1 / (1 + environ_var), digits = 2)
    )
  )
  dad.mean <- colSums(g$father_geno) / L
  mum.mean <- colSums(g$mother_geno) / L
  abline(lm(dad.mean ~ mum.mean), col = "red")
  text(x = 0.5,
       y = 1.5,
       format(cor(mum.mean, dad.mean), dig = 2),
       col = "red")



  ### Make a child

  ###CHECK THIS, WILL NEED TO DO TWICE ONCE FOR DAUGHTER ONCE FOR SONS?
  dad_gamete <- apply(g$father_geno, 2, make_a_gamete)
  mum_gamete <- apply(g$mother_geno, 2, make_a_gamete)
  child_genotype <- dad_gamete + mum_gamete
  cat(mean(child_genotype == 0),
      mean(child_genotype == 1),
      mean(child_genotype == 2),
      "\n")

  childs_phenotype_A <-
    generate_phenotype(eff = g$effects_A, child_genotype)
  childs_phenotype_B <-
    generate_phenotype(eff = g$effect_B, child_genotype)

  tmp1 <- (apply(child_genotype[loci_A, ], 1, function(test.geno) {
    summary(lm(childs_phenotype_A[, 1] ~ test.geno))$coeff[2, ]
  }))
  tmp2 <- (apply(child_genotype[loci_A,], 1, function(test.geno) {
    summary(lm(childs_phenotype_A[, 1] ~ test.geno))$coeff[2, ]
  }))
  tmp3 <- (apply(child_genotype[loci_B, ], 1, function(test.geno) {
    summary(lm(childs_phenotype_B[, 1] ~ test.geno))$coeff[2, ]
  }))
  tmp4 <- (apply(child_genotype[loci_B,], 1, function(test.geno) {
    summary(lm(childs_phenotype_B[, 1] ~ test.geno))$coeff[2, ]
  }))
  #	col <- c(rep("blue", length(loci_A)), rep("green", length(loci_B)))
  plot(
    tmp1[1, ],
    tmp2[1, ],
    xlab = "effect size trait 1",
    ylab = "effect size trait 2",
    cex.lab = 1.5
    ,
    pch = 19,
    col = "green",
    xlim = range(c(tmp1[1, ], tmp2[1, ])),
    ylim = range(c(tmp2[1, ], tmp4[1, ]))
  )
  points(
    tmp3[1, ],
    tmp4[1, ],
    xlab = "effect size trait 1",
    ylab = "effect size trait 2",
    cex.lab = 1.5
    ,
    pch = 19,
    col = "blue"
  )

  print(summary(lm(tmp1[1, ] ~ tmp2[1, ])))
  print(summary(lm(tmp3[1, ] ~ tmp4[1, ])))


}

#' Provide a quick report on the simulation of parental pairs and a child.
#'
#' This report is only for a single generation
#' @param am_corr_coeff The associative mating correlation coefficent
#' @param ntrait_loci  The number of loci for each trait
#' @param n_loci The number of loci to simulate.  
#' @param environ_var environmental variance.  A vector of length 1 or 2.
#' @param pop_size The Size of the populations to simulate.  Equal numbers of males and
#' females of this size.
#' @param generations The number of generations to run the simulation.
#' @param ntraining The size of the training data set
#' @param af_r_func A function of one variable (n) that returns a vector of allele frequencies of length n
#' @export
simulate_gen_assoc <- function(
  am_corr_coeff, 
  ntrait_loci,
  nloci,
  envir_var,
  pop_size,
  generations,
  ntraining,
  af_r_func = runif) {
  
  allele_freq <- af_r_func(nloci)                         
  if (length(envir_var) ==1) envir_var <- rep(envir_var, 2)
  ## simulate training dataset for lasso 
  training <- sample_initial_genotypes(ntraining, allele_freq)
  
  trait_loci <- sample(nloci, 2*ntrait_loci)
  trait_loci_A <- trait_loci[1:ntrait_loci]
  trait_loci_B <- trait_loci[(ntrait_loci+1):(2*ntrait_loci)]
  
  
  # get the  effects for the traits
  effects_A <- generate_effects(trait_loci_A, allele_freq, env_var = envir_var[1])
  effects_B <- generate_effects(trait_loci_B, allele_freq, env_var = envir_var[2])
  # and phenotypes for training
  training_pheno_A <- generate_phenotype(effects_A, training)$pheno
  training_pheno_B <- generate_phenotype(effects_B, training)$pheno
  # this is an example estimated by doing a bigger analysis that takes a long time
  fitted_lambda <- 0.02310998
  lasso_model_A <-  glmnet::glmnet(t(training), training_pheno_A, alpha = 1, lambda = fitted_lambda)
  lasso_model_B <-  glmnet::glmnet(t(training), training_pheno_B, alpha = 1, lambda = fitted_lambda)
  
  rm(training, training_pheno_A, training_pheno_B)
  
  ## A function to collect standard deviations and correlations
  # This one is more complicated as it collects all standard deviation and all correlations named
  summarise_phenotypes3 <- function(gg, effects_A, effects_B, lasso_model_A, lasso_model_B) {
    ga <- cbind(gg$male_genotype, gg$female_genotype)
    A <- setNames(generate_phenotype(effects_A, ga), c("phenoA", "addgenA"))
    lasso_A <- predict(lasso_model_A, t(ga))
    B  <- setNames(generate_phenotype(effects_B, ga), c("phenoB", "addgenB"))
    lasso_B <- predict(lasso_model_B, t(ga))
    cr <- cor(cbind(A, lassoA=lasso_A[, 1], B, lassoB=lasso_B[,1]))
    v <- var(cbind(A, lassoA=lasso_A[, 1], B, lassoB=lasso_B[,1]))
    
    nn <- outer(colnames(v), colnames(v), paste, sep="_")
    crr <- cr[upper.tri(cr, diag=FALSE)]
    names(crr) <- nn[upper.tri(cr, diag=FALSE)]
    c(sqrt(diag(v)), crr[c(4, 8, 13)])
  
  }
  
  ## ----simulate_generations---------------------
  
  go <- list(male_genotype = sample_initial_genotypes(pop_size, allele_freq),
             female_genotype = sample_initial_genotypes(pop_size, allele_freq))
  
  p <- list()
  p[[1]] <- summarise_phenotypes3(go, effects_A, effects_B, lasso_model_A, lasso_model_B)
  for (i in 1:generations) {
  #  if (i%%5 == 0) cat("generation", i, "\n")
    gn <- generation(go, effects_A, effects_B, am_corr_coeff)
    go <- gn
    p[[i+1]] <- summarise_phenotypes3(gn, effects_A, effects_B, lasso_model_A, lasso_model_B)
  }
  return(
    setNames(
      data.frame(matrix(unlist(p), nrow = 1+generations, byrow=TRUE)),
      names(p[[1]])
    )
  )
}



#' Provide a quick report on the simulation of parental pairs and a child.
#'
#' This report is only for a single generation
#' @param am_corr_coeff The associative mating correlation coefficent
#' @param ntrait_loci  The number of loci for each trait
#' @param n_loci The number of loci to simulate.  
#' @param environ_var environmental variance.  A vector of length 1 or 2.
#' @param pop_size The Size of the populations to simulate.  Equal numbers of males and
#' females of this size.
#' @param generations1 The number of generations to run the simulation before the break.
#' @param generations2 The number of generations after the two new polygenic scores are calculated
#' @param ntraining The size of the training data set
#' @param af_r_func A function of one variable (n) that returns a vector of allele frequencies of length n
#' @export
simulate_gen_assoc_recalc_PGS <- function(
  am_corr_coeff, 
  ntrait_loci,
  nloci,
  envir_var,
  pop_size,
  generations1,
  generations2,
  ntraining,
  af_r_func = runif) {
  ## A function to collect standard deviations and correlations
  # This one is more complicated as it collects all standard deviation and all correlations named
  summarise_phenotypes5 <- function(gg, effects_A, effects_B, lasso_model_A, lasso_model_B, lasso_model_A2=NULL, lasso_model_B2=NULL) {
    ga <- cbind(gg$male_genotype, gg$female_genotype)
    A <- setNames(generate_phenotype(effects_A, ga), c("phenoA", "addgenA"))
    lasso_A <- predict(lasso_model_A, t(ga))[, 1]
    B  <- setNames(generate_phenotype(effects_B, ga), c("phenoB", "addgenB"))
    lasso_B <- predict(lasso_model_B, t(ga))[, 1]
    
    if (!is.null(lasso_model_A2))
      lasso_A2 <- predict(lasso_model_A2, t(ga))[, 1]
    else
      lasso_A2 <- rep(NA, length(lasso_A))
    
    if (!is.null(lasso_model_B2))
      lasso_B2 <- predict(lasso_model_B2, t(ga))[, 1]
    else
      lasso_B2 <- rep(NA, length(lasso_A))
    
    xxx <- cbind(A, lassoA=lasso_A, B, lassoB=lasso_B, lassoA2 = lasso_A2, lassoB2=lasso_B2)
    cr <- cor(xxx)
    v <-  var(xxx)
    
    nn <- outer(colnames(v), colnames(v), paste, sep="_")
    crr <- cr[upper.tri(cr, diag=FALSE)]
    names(crr) <- nn[upper.tri(cr, diag=FALSE)]
    c(sqrt(diag(v)), crr)
  }
  
  allele_freq <- af_r_func(nloci)                         
  if (length(envir_var) ==1) envir_var <- rep(envir_var, 2)
  ## simulate training dataset for lasso 
  training <- sample_initial_genotypes(ntraining, allele_freq)
  
  trait_loci <- sample(nloci, 2*ntrait_loci)
  trait_loci_A <- trait_loci[1:ntrait_loci]
  trait_loci_B <- trait_loci[(ntrait_loci+1):(2*ntrait_loci)]
  
  # get the  effects for the traits
  effects_A <- generate_effects(trait_loci_A, allele_freq, env_var = envir_var[1])
  effects_B <- generate_effects(trait_loci_B, allele_freq, env_var = envir_var[2])
  # and phenotypes for training
  training_pheno_A <- generate_phenotype(effects_A, training)$pheno
  training_pheno_B <- generate_phenotype(effects_B, training)$pheno
  # this is an example estimated by doing a bigger analysis that takes a long time
  fitted_lambda <- 0.02310998
  lasso_model_A <-  glmnet::glmnet(t(training), training_pheno_A, alpha = 1, lambda = fitted_lambda)
  lasso_model_B <-  glmnet::glmnet(t(training), training_pheno_B, alpha = 1, lambda = fitted_lambda)
  
  rm(training, training_pheno_A, training_pheno_B)
  

  ## ----simulate_generations---------------------
  
  go <- list(male_genotype = sample_initial_genotypes(pop_size, allele_freq),
             female_genotype = sample_initial_genotypes(pop_size, allele_freq))
  
  p <- list()
  p[[1]] <- summarise_phenotypes5(go, effects_A, effects_B, lasso_model_A, lasso_model_B, NULL, NULL)
  for (i in 1:generations1) {
    #  if (i%%5 == 0) cat("generation", i, "\n")
    gn <- generation(go, effects_A, effects_B, am_corr_coeff)
    go <- gn
    p[[i+1]] <- summarise_phenotypes5(gn, effects_A, effects_B, lasso_model_A, lasso_model_B)
  }
  
  training2 <- cbind(gn$male_genotype, gn$female_genotype)
  training_pheno_A2 <- generate_phenotype(effects_A, training2)$pheno
  training_pheno_B2 <- generate_phenotype(effects_B, training2)$pheno
  
  lasso_model_A2 <-  glmnet::glmnet(t(training2), training_pheno_A2, alpha = 1, lambda = fitted_lambda)
  lasso_model_B2 <-  glmnet::glmnet(t(training2), training_pheno_B2, alpha = 1, lambda = fitted_lambda)
  
  rm(training2, training_pheno_B2, training_pheno_A2)
  
  for (i in (generations1+1):(generations1+generations2)) {
    gn <- generation(go, effects_A, effects_B, am_corr_coeff)
    go <- gn
    p[[i+1]] <- summarise_phenotypes5(gn, effects_A, effects_B, lasso_model_A, lasso_model_B, lasso_model_A2, lasso_model_B2)
  }
  return(
    setNames(
      data.frame(matrix(unlist(p), nrow = 1+generations1+generations2, byrow=TRUE)),
      names(p[[1]])
    )[,c(1:8, 12, 16, 21, 26, 35,36)]
  )
}
