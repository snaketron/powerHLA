data {
  int <lower=0> K; //alleles (categories)
  int <lower=0> B; //number of samples in y_rng
  int rng_gamma [B, K]; //samples
  int gamma [K]; // background HLA allele counts
  vector [K] alpha; // alpha priors
}

generated quantities {
  simplex [K] p_rng [B];
  simplex [K] p;
  vector [K] d_diff [B];
  vector [K] d_or [B];
  vector [K] d_lor [B];

  // estimate p of background population (a=0)
  p = dirichlet_rng(to_vector(gamma));

  // estimate p of each simulated sample
  for(b in 1:B) {
    p_rng[b] = dirichlet_rng(to_vector(rng_gamma[b,]) + alpha);
    d_or[b] = p ./ p_rng[b];
    d_lor[b] = log(d_or[b]);
    d_diff[b] = p - p_rng[b];
  }
}
