data {
  int<lower=0> N;                           // number of observations
  int<lower=0> K;                           // number of assays
  int<lower=1, upper=K> assay[N];           // assay index vector
  vector[N] timeVec;                        // predictor
  vector<lower=0>[N] nPositive;             // Number of positive samples
  vector<lower=0>[N] nNegative;             // Number of samples tested
}

parameters {
  // 
  real timeSlope;                     // mean slope
  real intercept;                     // mean intercept
  real<lower=0> slopeSigma;           // sd of the intercept
  real<lower=0> interceptSigma;       // sd of the intercept
  vector[K] assayIntercept;           // intercept of each assay 
  vector[K] assaySlope;               // intercept change of each assay 
}

transformed parameters{
  vector<lower=0, upper=1>[N] sensitivity;  // variable to contain % outcome
  sensitivity = inv_logit(assayIntercept[assay] + assaySlope[assay] .* timeVec);  // logistic regression model 
}

model {
  timeSlope ~ normal(0, 100);
  intercept ~ normal(0, 100);
  slopeSigma ~ gamma(4, 8);
  interceptSigma ~ gamma(4, 4);
  assaySlope ~ normal(timeSlope, slopeSigma);
  assayIntercept ~ normal(intercept, interceptSigma);
  sensitivity ~ beta(nPositive+1, nNegative+1);
}

