## THIS SCRIPT CONTAINS ALL THE STAN MODELS ###

Two_HT <- "data {
int<lower=1> I;               // # questions
int<lower=1> J;               // # persons
int<lower=1> N;               // # observations
int<lower=1, upper=I> ii[N];  // question for n
int<lower=1, upper=J> jj[N];  // person for n
int<lower=0, upper=1> y[N];   // correctness for n
vector[N] z;   // foil or target
}

parameters {
vector<lower=0, upper=1>[J] theta;
vector<lower=0, upper=1>[J] gamma;
real<lower=0> alpha_gamma;
real<lower=0> beta_gamma;
}

model {
vector[N] p;
vector[N] g;
theta ~ normal(0.5,1);         // informative true prior
gamma ~ beta(alpha_gamma, beta_gamma);

p <- theta[jj];
g <- fabs(z-gamma[jj]);
y ~ bernoulli(p + (1-p).*g);
}
"

Rasch <- "data {
  int<lower=1> I;               // # questions
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // question for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n

}
parameters {
  vector[I] alpha_beta;
  vector[J] theta;
  real mu_beta;
  real<lower=0> sigma_beta;
}

model {
vector[N] eta;
theta ~ normal(0,1);
alpha_beta ~ normal(0,1);
mu_beta ~ cauchy(0,5);
sigma_beta ~ cauchy(0,5);
eta <- theta[jj] - (mu_beta+sigma_beta*alpha_beta[ii]);
y ~ bernoulli_logit(eta);
}
"



GCMR <-
  "data {
int<lower=1> J;              // number of students
int<lower=1> I;              // number of questions
int<lower=1> N;              // number of observations
int<lower=1,upper=J> jj[N];  // student for observation n
int<lower=1,upper=I> ii[N];  // question for observation n
int<lower=0,upper=1> y[N];   // correctness for observation n
vector<lower=-1,upper=1>[N] z;   // foil or target
//real gamma_m; //mean of gamma prior
//real gamma_sd; //sd of gamma prior
//real gamma_a; //skewness of gamma prior
}

parameters {
vector[J] theta;      // ability of student j - mean ability
vector[I] alpha_beta;       // difficulty of question k
//vector[J] alpha_gamma; //gamma parameter
vector<lower=0,upper=1>[J] gamma;
// hyperparameters
real<lower=0> sigma_beta;
//real<lower=0> sigma_theta;
//real<lower=0> sigma_gamma;
real mu_beta;
real<lower=0> alpha_gamma;
real<lower=0> beta_gamma;

}

model {
vector[N] p;
vector[N] g;
theta ~ normal(0,1);         // informative true prior
alpha_beta ~ normal(0,1);
mu_beta ~ cauchy(0,5);
sigma_beta ~ cauchy(0,5);
gamma ~ beta(alpha_gamma, beta_gamma);

p <- inv_logit(theta[jj]-(mu_beta+sigma_beta*alpha_beta[ii]));
g <- fabs(z-gamma[jj]);
y ~ bernoulli(p + (1-p).*g);
}
"



GCMR_Extended <-
  "data {
int<lower=1> J;              // number of students
int<lower=1> I;              // number of questions
int<lower=1> N;              // number of observations
int<lower=1,upper=J> jj[N];  // student for observation n
int<lower=1,upper=I> ii[N];  // question for observation n
int<lower=0,upper=1> y[N];   // correctness for observation n
real age[N];        // subject age for observation n
real<lower=-1,upper=1> domain[N];     // item domain (object vs. scene) for observation n
int unique_ii[I];
real unique_domain[I];
int unique_jj[J];
real unique_age[J];
vector<lower=-1,upper=1>[N] z;   // foil or target
}

parameters {
vector[J] theta; // ability of student j - mean ability
vector[I] alpha_beta; // difficulty of question k
vector<lower=0,upper=1>[J] gamma;
// hyperparameters
real<lower=0> sigma_beta;
real mu_beta;
real<lower=0> alpha_gamma;
real<lower=0> beta_gamma;
real a; // age regression coeff
real d; // domain regression coeff
real age_domain; //age by domain interaction coeff
}


model {
vector[N] p;
vector[N] g;
mu_beta ~ cauchy(0,5);
sigma_beta ~ cauchy(0,5);
gamma ~ beta(alpha_gamma, beta_gamma);
a ~ normal(0,0.215);
d ~ normal(0,0.504);
age_domain ~ normal(0,0.43);

for (i in 1:I){
alpha_beta[unique_ii[i]] ~ normal(d*unique_domain[i],1);
}

for (j in 1:J){
theta[unique_jj[j]] ~ normal(a*unique_age[j],1); 
}

for (n in 1:N){
p[n] = inv_logit(theta[jj[n]]-(mu_beta+sigma_beta*alpha_beta[ii[n]])+age_domain*domain[n]*age[n]);
g[n] = fabs(z[n]-gamma[jj[n]]);
y[n] ~ bernoulli(p[n] + [(1-p[n])]*g[n]);
}
}
"


