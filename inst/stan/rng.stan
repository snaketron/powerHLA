data {
  int<lower=0> K; // number of categories
  int<lower=0> N; // population (number of alleles) size
  simplex [K] af; // simplex of allele frequencies
}

generated quantities {
  int y [K];
  y = multinomial_rng(af, N); // generate allele counts
}
