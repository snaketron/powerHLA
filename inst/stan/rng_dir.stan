data {
  int<lower=0> K; // number of categories
  vector <lower=0> [K] alpha; // alphas
}

generated quantities {
  simplex [K] theta;
  theta = dirichlet_rng(alpha); // generate allele frequencies
}
