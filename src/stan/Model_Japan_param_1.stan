Model_1 = "
functions{
vector convolution(vector X, vector Yrev, int K) {
vector[K-1] res;
res[1] = X[1]*Yrev[K];
for (k in 2:K-1) // 2:K-1 is equivalent to 2:(K-1)
res[k] = dot_product(head(X, k), tail(Yrev, k)); // by definition of the convolution
   return res;        
    }
}

data {
// number of week
int T;
// population
int N;
int num_data;
// delay
int l;
int delay;
// unvaccinated case
vector[T+l] it;
// vaccinated case
vector[T+l] jt;
// vaccination
vector[T+l+delay] Gamma;
// serial interval
real SI[T+l];
// odds_
vector[T+l+delay] odds;
real others[T+l+delay];
//VOC alpha;
real alpha[T+l+delay];
//VOC delta;
real delta[T+l+delay];
//VOC delta;
real omicron[T+l+delay];
}

transformed data{
vector[T+l] SI_rev;
vector[T+l-1] conv;
for(h in 1:T+l)
SI_rev[h] = SI[T+l-h+1];

conv = convolution(it+jt, SI_rev, T+l);
}

parameters{
real<lower=0> k[4];
vector<lower=0>[T-1] Rit;
vector<lower=0,upper=1>[T-1] eps;
real<lower=0> eta[1];
real<lower=0,upper=1> c[4];
}

transformed parameters{
real<lower=0,upper=1> zeta[T+l+delay-1];
vector<lower=0>[T-1] Rjt;

for(s in 1:T+l+delay-1){
real ve_reduction_o[s];
real ve_reduction_a[s];
real ve_reduction_d[s];
real ve_reduction_om[s];
real vax_rev[s];
real convolution_r[s];

for(t in 1:s){                                   
ve_reduction_o[t] = c[1] * exp(-k[1]*(t-1));
ve_reduction_a[t] = c[2] * exp(-k[2]*(t-1));
ve_reduction_d[t] = c[3] * exp(-k[3]*(t-1));
ve_reduction_om[t] = c[4] * exp(-k[4]*(t-1));
vax_rev[t] = Gamma[s-t+1]; 
convolution_r[t] = (others[s] * ve_reduction_o[t] + alpha[s] * ve_reduction_a[t] + delta[s] * ve_reduction_d[t] + omicron[s] * ve_reduction_om[t]) * vax_rev[t];
}
zeta[s] = sum(convolution_r);
}
for(t in 1:T-1)
Rjt[t] = odds[t+l+delay]  * (1-eps[t]) * Rit[t];
}

model{ 
Rit ~ normal(0.5,1);
k ~ normal(0,10);
c ~ beta(5,2);
eta ~ normal(0,100);
for(t in 1:T-1)
eps[t] ~ beta((eta[1]*jt[t+l+1]) *zeta[t+l+delay],(eta[1]*jt[t+l+1]) -(eta[1]*jt[t+l+1]) *zeta[t+l+delay]);
target += gamma_lpdf(it[1+l+1:T+l] | Rit .* conv[1+l:T+l-1] + 1e-13, 1.0) + gamma_lpdf(jt[1+l+1:T+l] | Rjt .* conv[1+l:T+l-1] + 1e-13, 1.0);
}

generated quantities{
real ve_o[num_data];
real ve_a[num_data];
real ve_d[num_data];
real ve_om[num_data];
real ii[T-1];
real jj[T-1];

for(t in 1:num_data){
ve_o[t] = c[1] * exp(-k[1]*(t-1));
ve_a[t] = c[2] * exp(-k[2]*(t-1));
ve_d[t] = c[3] * exp(-k[3]*(t-1));
ve_om[t] = c[4] * exp(-k[4]*(t-1));
}
for(t in 1:T-1){
ii[t] = Rit[t] * conv[t+l];
jj[t] = Rjt[t] * conv[t+l];
}
}
"