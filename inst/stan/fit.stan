data {
  int<lower=0> N;
  vector<lower=0> [N] alpha;
  int<lower=0> gamma [N];
}

parameters {
  simplex [N] theta;
}

model {
  gamma ~ multinomial(theta);
  theta ~ dirichlet(alpha);
}

generated quantities {
  vector [N] theta_rng = dirichlet_rng(to_vector(gamma) + alpha);
}
