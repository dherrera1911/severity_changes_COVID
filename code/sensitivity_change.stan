data {
  int<lower=0> N;                           // number of observations
  int<lower=0> K;                           // number of assays
  int<lower=0> M;                           // number of studies:assays
  int<lower=1, upper=K> assay[N];           // assay index vector
  int<lower=1, upper=M> study[N];           // study:assay index vector
  vector[N] timeVec;                        // predictor
  int<lower=0> nPositive[N];             // Number of positive samples
  int<lower=0> nTested[N];             // Number of samples tested
}

parameters {
  // 
  real timeSlope;                     // mean slope
  real intercept;                     // mean intercept
  real<lower=0> slopeSigma;           // sd of the intercept
  real<lower=0> interceptSigma;       // sd of the intercept
  real<lower=0> studySigma;           // sd of the study effect
  vector[K] assayIntercept;           // intercept of each assay 
  vector[K] assaySlope;               // intercept change of each assay 
  vector[M] studyIntercept;           // intercept of each study:assay
}

transformed parameters{
  vector<lower=0, upper=1>[N] sensitivity;  // variable to contain % outcome
  sensitivity = inv_logit(assayIntercept[assay] + studyIntercept[study] + assaySlope[assay] .* timeVec);  // logistic regression model 
}

model {
  timeSlope ~ normal(0, 100);
  intercept ~ normal(0, 100);
  slopeSigma ~ gamma(4, 8);
  interceptSigma ~ gamma(4, 4);
  studySigma ~ gamma(4, 4);
  assaySlope ~ normal(timeSlope, slopeSigma);
  assayIntercept ~ normal(intercept, interceptSigma);
  studyIntercept ~ normal(0, studySigma);
  nPositive ~ binomial(nTested, sensitivity);
}


