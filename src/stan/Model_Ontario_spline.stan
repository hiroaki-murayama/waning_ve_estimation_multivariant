Model_s = "
functions{
vector convolution(vector X, vector Yrev, int K) {
vector[K-1] res;
res[1] = X[1]*Yrev[K];
for (k in 2:K-1) // 2:K-1 is equivalent to 2:(K-1)
res[k] = dot_product(head(X, k), tail(Yrev, k)); // by definition of the convolution
return res;        
}

// get the vector of spacings between nodes
  vector spline_geths(int n_nodes, vector nodes){
int n = n_nodes -1;
vector[n] hs;
for (i in 1:n){
hs[i] = nodes[i+1]-nodes[i];
}
return hs;
}
  
// obtain the vector of spline coefficients given the location
// of the nodes and values there
// We are using natural spline definition           
vector spline_getcoeffs(int n_nodes, vector nodes, vector vals){
int n=n_nodes-1;
vector[n] hi;
vector[n] bi;
vector[n-1] vi;
vector[n-1] ui;
vector[n_nodes] ret;
vector[n-1] zs;
matrix[n-1,n-1] M = rep_matrix(0, n-1, n-1);

n = n_nodes-1;

for (i in 1:n){
hi[i] = nodes[i+1]-nodes[i];
bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
}
for (i in 2:n){
vi[i-1] = 2*(hi[i-1]+hi[i]);
ui[i-1] = 6*(bi[i] - bi[i-1]);
}
for (i in 1:n-1)
M[i,i] = vi[i];
for (i in 1:n-2){
M[i+1,i] = hi[i+1];
M[i,i+1] = hi[i+1];
}
zs = inverse(M) * ui;
    
ret[1]=0;
ret[n_nodes] =0;
ret[2:n_nodes-1]=zs;

return ret;
}

// Evaluate the spline, given nodes, values at the nodes
// spline coefficients, locations of evaluation points
// and integer bin ids of each point            
vector spline_eval(int n_nodes, vector nodes, vector vals, vector zs, int n_dat, vector x, int[] i){

vector[n_nodes-1] h;
vector[n_dat] ret;
int i1[n_dat];
for (ii in 1:n_dat){
i1[ii] = i[ii] + 1;
}
h = spline_geths(n_nodes, nodes);

    ret = (
           zs[i1] ./ 6 ./ h[i] .* square(x-nodes[i]) .* (x-nodes[i])+
           zs[i]  ./ 6 ./ h[i] .* square(nodes[i1]-x) .* (nodes[i1]-x)+
           (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6) .* (x-nodes[i])+
           (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6) .* (nodes[i1]-x)
           );
return ret;
}

// find in which node interval we should place each point of the vector                   
int[] spline_findpos(int n_nodes, vector nodes, int n_dat, vector x){
int ret[n_dat];
for (i in 1:n_dat){
for (j in 1:n_nodes-1){
if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1])){
ret[i] = j;
}
if (x[i]==nodes[n_nodes]){
ret[i] = j;
}
}
}
return ret;
}
}


data {
// number of week
int T;
// population in New York City
int N;
int l;
// delay
int delay;
// unvaccinated case
vector[T+l] it;
// vaccinated case
vector[T+l] jt;
// vaccination
vector[T+l+delay] Gamma;
// serial interval
real SI[T+l];
// odds
vector[T+l+delay] odds;
//VOC others;
real others[T+l+delay];
//VOC alpha;
real alpha[T+l+delay];
//VOC delta;
real delta[T+l+delay];
int num_data;   // number of data points
int nknots;
vector[num_data] X;
vector[nknots] xknots;
}

transformed data{
vector[T+l] SI_rev;
vector[T+l-1] conv;
// determine which knots the point belong to
int x_pos_knots[num_data] = spline_findpos(nknots, xknots, num_data, X);

for(h in 1:T+l)
SI_rev[h] = SI[T+l-h+1];

conv = convolution(it+jt, SI_rev, T+l);
}

parameters{
vector<lower=0>[T-1] Rit;
// the parameters of our spline model are the values at the knots
//vector[nknots] yknots_om;
vector[nknots] yknots_o;
vector[nknots] yknots_a;
vector[nknots] yknots_d;
vector<lower=0,upper=1>[T-1] eps;
real<lower=0> eta[1];
}

transformed parameters{
real zeta[T+l+delay-1];
vector<lower=0>[T-1] Rjt;
//vector[nknots] spl_coeffs_om = spline_getcoeffs(nknots, xknots, yknots_om);
vector[nknots] spl_coeffs_o = spline_getcoeffs(nknots, xknots, yknots_o);
vector[nknots] spl_coeffs_a = spline_getcoeffs(nknots, xknots, yknots_a);
vector[nknots] spl_coeffs_d = spline_getcoeffs(nknots, xknots, yknots_d);
//vector[num_data] nc_om = spline_eval(nknots, xknots, yknots_om, spl_coeffs_om, num_data, X, x_pos_knots);
vector[num_data] nc_o = spline_eval(nknots, xknots, yknots_o, spl_coeffs_o, num_data, X, x_pos_knots);
vector[num_data] nc_a = spline_eval(nknots, xknots, yknots_a, spl_coeffs_a, num_data, X, x_pos_knots);
vector[num_data] nc_d = spline_eval(nknots, xknots, yknots_d, spl_coeffs_d, num_data, X, x_pos_knots);

for(s in 1:(T+l+delay-1)){
//real ve_reduction_om[s];
real ve_reduction_o[s];
real ve_reduction_a[s];
real ve_reduction_d[s];
real vax_rev[s];
real convolution_r[s];
for(t in 1:s){
ve_reduction_o[t] = inv_logit(nc_o[t]);
ve_reduction_a[t] = inv_logit(nc_a[t]);
ve_reduction_d[t] = inv_logit(nc_d[t]);
vax_rev[t] = Gamma[s-t+1]; 
convolution_r[t] = (others[s] * ve_reduction_o[t] + alpha[s] * ve_reduction_a[t] + delta[s] * ve_reduction_d[t]) * vax_rev[t];
}
zeta[s] = sum(convolution_r);
}

for(t in 1:T-1)
Rjt[t] =  odds[t+l+delay] * (1-eps[t]) * Rit[t];
}

model{ 
for(t in 1:T-1){
eps[t] ~ beta((eta[1]/sqrt(jt[t+l])) *zeta[t+l+delay],(eta[1]/sqrt(jt[t+l]))-(eta[1]/sqrt(jt[t+l]))*zeta[t+l+delay]);
}
target += gamma_lpdf(it[1+l:T+l-1] | Rit .* conv[1+l:T+l-1] + 1e-13, 1.0) + gamma_lpdf(jt[1+l:T+l-1] | Rjt .* conv[1+l:T+l-1] + 1e-13, 1.0);

inv_logit(yknots_o) ~ beta(5,2);
inv_logit(yknots_a) ~ beta(5,2);
inv_logit(yknots_d) ~ beta(5,2);
Rit ~ normal(0.5,1);
eta ~ normal(0,200);
}

generated quantities{
real ve_o[num_data];
real ve_a[num_data];
real ve_d[num_data];
real ii[T-1];
real jj[T-1];
for(t in 1:num_data){
//ve_reduction_om[t] = inv_logit(nc_om[t]);
ve_o[t] = inv_logit(nc_o[t]);
ve_a[t] = inv_logit(nc_a[t]);
ve_d[t] = inv_logit(nc_d[t]);

}
for(t in 1:T-1){
ii[t] = Rit[t] * conv[t+l];
jj[t] = Rjt[t] * conv[t+l];
}
}
"