data {
  int<lower=0> K; // number of categories
  int<lower=0> N; // population (number of alleles) size
  simplex [K] theta; // simplex of allele frequencies
}

generated quantities {
  int gamma [K];
  gamma = multinomial_rng(theta, N); // generate allele counts
}
