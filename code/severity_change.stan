functions{
  real bin_lpmf(int x, real n, real theta){
    real ans  = lchoose(n, x) + x*log(theta) + (n-x)*log1m(theta);
    return(ans);
  }
}

data {
  int<lower=0> N;                       // number of observations
  int<lower=0> K;                       // number of locations
  int<lower=1, upper=K> location[N];    // location index vector
  vector[N] ageVec;                     // predictor
  vector<lower=0>[N] population;                // population count
  vector<lower=0>[N] seroprev1Shape1;           // Prevalence beta shape1
  vector<lower=0>[N] seroprev1Shape2;           // Prevalence beta shape2
  vector<lower=0>[N] seroprev2Shape1;           // Prevalence beta shape1
  vector<lower=0>[N] seroprev2Shape2;           // Prevalence beta shape2
  int<lower=0> outcomes1[N];                    // outcome count
  int<lower=0> outcomes2[N];                    // outcome count
}

parameters {
  // Countries outcome fit
  real ageSlope;                      // mean slope
  real intercept;                     // mean intercept
  real interceptChange;               // mean effect of period
  real<lower=0> interceptSigma;       // sd of the intercept
  real<lower=0> changeSigma;       // sd of the change in time
  vector[K] locationIntercept;        // intercept of each location 
  vector[K] locationChange;        // intercept change of each location 
  vector<lower=0, upper=1>[N] prevalence_raw1; 
  vector<lower=0, upper=1>[N] prevalence_diff; 
}

transformed parameters{
  vector<lower=0, upper=1>[N] outcomeRate1;  // variable to contain % outcome
  vector<lower=0, upper=1>[N] outcomeRate2;  // variable to contain % outcome
  vector<lower=0, upper=1>[N] prevalence_raw2;
  outcomeRate1 = inv_logit(locationIntercept[location] + ageSlope * ageVec);  // logistic regression model 
  outcomeRate2 = inv_logit(locationIntercept[location] + locationChange[location] + ageSlope * ageVec);  // logistic regression model 
  prevalence_raw2 = prevalence_diff + prevalence_raw1;
}

model {
  ageSlope ~ normal(0, 100);
  intercept ~ normal(0, 100);
  interceptChange ~ normal(0, 100);
  interceptSigma ~ gamma(4, 4);
  changeSigma ~ gamma(4, 10);
  locationIntercept ~ normal(intercept, interceptSigma);
  locationChange ~ normal(interceptChange, changeSigma);
  prevalence_raw1 ~ beta(seroprev1Shape1, seroprev1Shape2);
  prevalence_raw2 ~ beta(seroprev2Shape1, seroprev2Shape2);
  for (n in 1:N) {
    outcomes1[n] ~ bin(population[n] * prevalence_raw1[n], outcomeRate1[n]);
    outcomes2[n] ~ bin(population[n] * prevalence_diff[n], outcomeRate2[n]);
  }
}

