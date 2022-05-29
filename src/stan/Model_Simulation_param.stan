Model = "
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
// population in New York city
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
real Gamma[T+l+delay];
// serial interval
real SI[T+l];
// odds
vector[T+l+delay] odds;
//Variant A;
real vA[T+l+delay];
//Variant B;
real vB[T+l+delay];
//Variant C;
real vC[T+l+delay];
//Variant D ;
real vD[T+l+delay];
}

transformed data{
vector[T+l] SI_rev;
vector[T+l-1] conv;
for(h in 1:T+l)
SI_rev[h] = SI[T+l-h+1];

conv = convolution(it+jt, SI_rev, T+l);
}

parameters{
real<lower=1,upper=2> p[4];
real<lower=0> k[4];
vector<lower=0>[T-1] Rit;
vector<lower=0,upper=1>[T-1] eps;
real<lower=0> eta[1];
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
ve_reduction_o[t] = p[1] * (1-inv_logit(k[1]*(t-1)));
ve_reduction_a[t] = p[2] * (1-inv_logit(k[2]*(t-1)));
ve_reduction_d[t] = p[3] * (1-inv_logit(k[3]*(t-1)));
ve_reduction_om[t] = p[4] * (1-inv_logit(k[4]*(t-1)));
vax_rev[t] = Gamma[s-t+1]; 
convolution_r[t] = (vA[s] * ve_reduction_o[t] + vB[s] * ve_reduction_a[t] + vC[s] * ve_reduction_d[t] + vD[s] * ve_reduction_om[t]) * vax_rev[t];
}
zeta[s] = sum(convolution_r);
}
for(t in 1:T-1)
Rjt[t] = odds[t+l+delay] * (1-eps[t]) * Rit[t];
}

model{ 
for(t in 1:T-1)
eps[t] ~ beta((eta[1]/sqrt(jt[t+l]+1e-3)) *zeta[t+l+delay],(eta[1]/sqrt(jt[t+l]+1e-3))-(eta[1]/sqrt(jt[t+l]+1e-3))*zeta[t+l+delay]);

Rit ~ normal(0.5,1);
k ~ normal(0,10);
eta ~ normal(0,100);
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
ve_o[t] = p[1] * (1-inv_logit(k[1]*(t-1)));
ve_a[t] = p[2] * (1-inv_logit(k[2]*(t-1)));
ve_d[t] = p[3] * (1-inv_logit(k[3]*(t-1)));
ve_om[t] = p[4] * (1-inv_logit(k[4]*(t-1)));
}
for(t in 1:T-1){
ii[t] = Rit[t] * conv[t+l];
jj[t] = Rjt[t] * conv[t+l];
}
}
"