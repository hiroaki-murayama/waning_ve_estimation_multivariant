set.seed(12345)
T=nrow(ncov_vax);
l = 2; 
N = 15000000#14826276;
Ri_A = rep(0,T-l); Ri_B = rep(0,T-l); Ri_C = rep(0,T-l); Ri_D = rep(0,T-l);
Rj_A = rep(0,T-l); Rj_B = rep(0,T-l); Rj_C = rep(0,T-l); Rj_D = rep(0,T-l); 
i0_A <- c(5300,5200); j0_A <- c(30,50); 
i0_B = c(20); j0_B = c(40); i0_C = c(20); j0_C = c(50); i0_D = c(50); j0_D = c(100);
eta <- 0.2
zeta_A <- rep(0,T-l); zeta_B <- rep(0,T-l); zeta_C <- rep(0,T-l); zeta_D <- rep(0,T-l)
eps_A <- rep(0,T-l); eps_B <- rep(0,T-l); eps_C <- rep(0,T-l); eps_D <- rep(0,T-l)
rollout <- ncov_vax$cumulative
rollout1 <- ncov_vax$vaccine/N
se <- 0.1
GT <- function(t){plnorm(t*7, log(4.699), log(2.936))-plnorm(t*7-7, log(4.699), log(2.936))}   
it_A <- rep(0, T); it_B <- rep(0, T); it_C <- rep(0, T); it_D <- rep(0, T);
jt_A <- rep(0, T); jt_B <- rep(0, T); jt_C <- rep(0, T); jt_D <- rep(0, T);
pA <- c(1,1,rep(0, T-l-l)); pB <- rep(0, T-l); pC <- rep(0, T-l); pD <- rep(0, T-l);

iit_A <- c(i0_A,rep(0, T-l)); iit_B <- rep(0, T); iit_C <- rep(0, T); iit_D <- rep(0, T) 
jjt_A <- c(j0_A,rep(0, T-l)); jjt_B <- rep(0, T); jjt_C <- rep(0, T); jjt_D <- rep(0, T)

iit <- c(i0_A,rep(0, T-l)); jjt <- c(j0_A,rep(0, T-l));

NPIs=0.4; NPIs1=0.4*0.3


for(t in (1+l):T){
  fn_overall_A <- rep(0,t); fn_overall_B <- rep(0,t); fn_overall_C <- rep(0,t); fn_overall_D <- rep(0,t);
  for(tt in 1:t){
    fn_overall_A[tt] <- ve_A(tt) * rollout1[t-tt+1] 
    fn_overall_B[tt] <- ve_B(tt) * rollout1[t-tt+1] 
    fn_overall_C[tt] <- ve_C(tt) * rollout1[t-tt+1] 
    fn_overall_D[tt] <- ve_D(tt) * rollout1[t-tt+1] 
  }   
  pA[t-l] = (iit_A[t-1] + jjt_A[t-1])/(iit[t-1] + jjt[t-1]); pB[t-l] = (iit_B[t-1] + jjt_B[t-1])/(iit[t-1] + jjt[t-1]); pC[t-l] = (iit_C[t-1] + jjt_C[t-1])/(iit[t-1] + jjt[t-1]); pD[t-l] = (iit_D[t-1] + jjt_D[t-1])/(iit[t-1] + jjt[t-1])
  
  zeta_A[t-l] <- sum(fn_overall_A[1:t])
  zeta_B[t-l] <- sum(fn_overall_B[1:t])
  zeta_C[t-l] <- sum(fn_overall_C[1:t])
  zeta_D[t-l] <- sum(fn_overall_D[1:t])
  eps_A[t-l] <- rbeta(1, shape1 = zeta_A[t-l] * eta * (jjt_A[t-1]), shape2 = eta * (jjt_A[t-1]) - zeta_A[t-l] * eta * (jjt_A[t-1]))
  eps_B[t-l] <- rbeta(1, shape1 = zeta_B[t-l] * eta * (jjt_B[t-1]), shape2 = eta * (jjt_B[t-1]) - zeta_B[t-l] * eta * (jjt_B[t-1]))
  eps_C[t-l] <- rbeta(1, shape1 = zeta_C[t-l] * eta * (jjt_C[t-1]), shape2 = eta * (jjt_C[t-1]) - zeta_C[t-l] * eta * (jjt_C[t-1]))
  eps_D[t-l] <- rbeta(1, shape1 = zeta_D[t-l] * eta * (jjt_D[t-1]), shape2 = eta * (jjt_D[t-1]) - zeta_D[t-l] * eta * (jjt_D[t-1]))
  
  
  Ri_A[t-l] <- exp(rnorm(1,log(0.5),se))
  Ri_B[t-l] <- exp(rnorm(1,log(0.5*1.5),se))
  Ri_C[t-l] <- exp(rnorm(1,log(0.5*1.5*2),se))
  Ri_D[t-l] <- exp(rnorm(1,log(0.5*1.5*2*2),se))
  Rj_A[t-l] <- (1-eps_A[t-l]) * (rollout[t]/(N-rollout[t])) * Ri_A[t-l] 
  Rj_B[t-l] <- (1-eps_B[t-l]) * (rollout[t]/(N-rollout[t])) * Ri_B[t-l] 
  Rj_C[t-l] <- (1-eps_C[t-l]) * (rollout[t]/(N-rollout[t])) * Ri_C[t-l] 
  Rj_D[t-l] <- (1-eps_D[t-l]) * (rollout[t]/(N-rollout[t])) * Ri_D[t-l] 
  
  
  if(t>=(1+l) & t<=14){
    for (tau in 2:t){it_A[tau-1] = Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
  } 
  if(t==15){
    for (tau in 2:t){it_A[tau-1] = Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = i0_B; jjt_B[t] = j0_B
  } 
  if(t>=16 & t<=29){
    for (tau in 2:t){it_A[tau-1] = Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
  } 
  if(t==30){
    for (tau in 2:t){it_A[tau-1] = NPIs * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = NPIs * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = NPIs * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = NPIs * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = i0_C; jjt_C[t] = j0_C;
  } 
  if(t>=31 & t<=54){
    for (tau in 2:t){it_A[tau-1] = NPIs * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = NPIs * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = NPIs * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = NPIs * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = NPIs * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = NPIs * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
  } 
  if(t==55){
    for (tau in 2:t){it_A[tau-1] = NPIs1 * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = NPIs1 * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = NPIs1 * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = NPIs1 * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = NPIs1 * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = NPIs1 * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
    iit_D[t] = i0_D; jjt_D[t] = j0_D;
  } 
  if(t>=56 & t<=T){
    for (tau in 2:t){it_A[tau-1] = NPIs1 * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = NPIs1 * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = NPIs1 * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = NPIs1 * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = NPIs1 * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = NPIs1 * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_D[tau-1] = NPIs1 * Ri_D[t-l] * (iit_D[t-tau+1] + jjt_D[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_D[tau-1] = NPIs1 * Rj_D[t-l] * (iit_D[t-tau+1] + jjt_D[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]))
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
    iit_D[t] = rpois(1,sum(it_D[1:(t-1)])); jjt_D[t] = rpois(1,sum(jt_D[1:(t-1)]));
  }
  iit[t] = iit_A[t] + iit_B[t] + iit_C[t] + iit_D[t] ## reproductive property of Poisson distr
  jjt[t] = jjt_A[t] + jjt_B[t] + jjt_C[t] + jjt_D[t] # pA[is.na(pA)] <- 0
}
