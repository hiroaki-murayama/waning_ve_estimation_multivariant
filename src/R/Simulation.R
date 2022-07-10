#### simulation code
set.seed(12345)
T=nrow(ncov_vax);
l = 2; 
N = 15000000#14826276;
Ri_A = rep(0,T-l); Ri_B = rep(0,T-l); Ri_C = rep(0,T-l); Ri_D = rep(0,T-l);
Rj_A = rep(0,T-l); Rj_B = rep(0,T-l); Rj_C = rep(0,T-l); Rj_D = rep(0,T-l); 
i0_A <- c(3300,3200); j0_A <- c(3,5); 
i0_B = c(30); j0_B = c(15); i0_C = c(30); j0_C = c(30); i0_D = c(5); j0_D = c(15);
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

NPIs=0.5; NPIs1=0.5*0.3


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
  
  
  Ri_A[t-l] <- exp(rnorm(1,log(1.3),se))
  Ri_B[t-l] <- exp(rnorm(1,log(1.3*1.85),se))
  Ri_C[t-l] <- exp(rnorm(1,log(1.3*1.85*2.1),se))
  Ri_D[t-l] <- exp(rnorm(1,log(1.3*1.85*2.1*3),se))
  Rj_A[t-l] <- (1-eps_A[t-l]) * Ri_A[t-l] 
  Rj_B[t-l] <- (1-eps_B[t-l]) * Ri_B[t-l] 
  Rj_C[t-l] <- (1-eps_C[t-l]) * Ri_C[t-l] 
  Rj_D[t-l] <- (1-eps_D[t-l]) * Ri_D[t-l] 
  
  
  if(t>=(1+l) & t<=14){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
  } 
  if(t==15){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = i0_B; jjt_B[t] = j0_B
  } 
  if(t>=16 & t<=29){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = (1-rollout[t-1]/N) * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = (rollout[t-1]/N) * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
  } 
  if(t==30){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * NPIs * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * NPIs * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = (1-rollout[t-1]/N) * NPIs * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = (rollout[t-1]/N) * NPIs * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = i0_C; jjt_C[t] = j0_C;
  } 
  if(t>=31 & t<=54){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * NPIs * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * NPIs * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = (1-rollout[t-1]/N) * NPIs * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = (rollout[t-1]/N) * NPIs * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = (1-rollout[t-1]/N) * NPIs * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = (rollout[t-1]/N) * NPIs * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
  } 
  if(t==55){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]));
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
    iit_D[t] = i0_D; jjt_D[t] = j0_D;
  } 
  if(t>=56 & t<=T){
    for (tau in 2:t){it_A[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_A[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_A[t-l] * (iit_A[t-tau+1] + jjt_A[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_B[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_B[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_B[t-l] * (iit_B[t-tau+1] + jjt_B[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_C[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_C[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_C[t-l] * (iit_C[t-tau+1] + jjt_C[t-tau+1]) * GT(tau-1)} 
    for (tau in 2:t){it_D[tau-1] = (1-rollout[t-1]/N) * NPIs1 * Ri_D[t-l] * (iit_D[t-tau+1] + jjt_D[t-tau+1]) * GT(tau-1)}
    for (tau in 2:t){jt_D[tau-1] = (rollout[t-1]/N) * NPIs1 * Rj_D[t-l] * (iit_D[t-tau+1] + jjt_D[t-tau+1]) * GT(tau-1)} 
    iit_A[t] = rpois(1,sum(it_A[1:(t-1)])); jjt_A[t] = rpois(1,sum(jt_A[1:(t-1)]))
    iit_B[t] = rpois(1,sum(it_B[1:(t-1)])); jjt_B[t] = rpois(1,sum(jt_B[1:(t-1)])); 
    iit_C[t] = rpois(1,sum(it_C[1:(t-1)])); jjt_C[t] = rpois(1,sum(jt_C[1:(t-1)])); 
    iit_D[t] = rpois(1,sum(it_D[1:(t-1)])); jjt_D[t] = rpois(1,sum(jt_D[1:(t-1)]));
  }
  iit[t] = iit_A[t] + iit_B[t] + iit_C[t] + iit_D[t] ## reproductive property of Poisson distr
  jjt[t] = jjt_A[t] + jjt_B[t] + jjt_C[t] + jjt_D[t] # pA[is.na(pA)] <- 0
}


#### arrange the output of simulation

### dataframe of case by vaccination status and variant
df_unvac <- rbind(iit_A,iit_B,iit_C,iit_D) %>% t()
df_vac <- rbind(jjt_A,jjt_B,jjt_C,jjt_D) %>% t()
label <- c(rep("A",T),rep("B",T),rep("C",T),rep("D",T)) %>% as.data.frame()
week <- seq(1,T,1) %>% rep(4) %>% as.data.frame()
case_df_type <- cbind(week,df_unvac,df_vac,label)
colnames(case_df_type) <- c("Week","Unvaccinated_A","Unvaccinated_B","Unvaccinated_C","Unvaccinated_D","Vaccinated_A","Vaccinated_B","Vaccinated_C","Vaccinated_D","Variant")

### dataframe of case by variant for the analysis
A_prop <- c(1,(iit_A[1:(T-1)]+jjt_A[1:(T-1)])/(iit[1:(T-1)]+jjt[1:(T-1)])) %>% as.data.frame()
B_prop <- c(0,(iit_B[1:(T-1)]+jjt_B[1:(T-1)])/(iit[1:(T-1)]+jjt[1:(T-1)])) %>% as.data.frame()
C_prop <- c(0,(iit_C[1:(T-1)]+jjt_C[1:(T-1)])/(iit[1:(T-1)]+jjt[1:(T-1)])) %>% as.data.frame()
D_prop <- c(0,(iit_D[1:(T-1)]+jjt_D[1:(T-1)])/(iit[1:(T-1)]+jjt[1:(T-1)])) %>% as.data.frame()
df_voc <- rbind(A_prop, B_prop, C_prop, D_prop) 
label <- c(rep("Variant A",T),rep("Variant B",T),rep("Variant C",T),rep("Variant D",T)) %>% as.data.frame()
week <- seq(1,T,1) %>% rep(4) %>% as.data.frame()
case_voc <- cbind(week,df_voc,label)
colnames(case_voc) <- c("Week","Proportion","Variant")

A <- case_voc %>% filter(Variant == "Variant A")
B <- case_voc %>% filter(Variant == "Variant B")
C <- case_voc %>% filter(Variant == "Variant C")
D <- case_voc %>% filter(Variant == "Variant D")



### dataframe of case by vaccination status for the analysis
df_unvac_non <- (iit[1:T]) %>% as.data.frame()
df_vac_non <- (jjt[1:T]) %>% as.data.frame()
week1 <- seq(1,T,1)
case_df <- cbind(week1,df_unvac_non,df_vac_non)
colnames(case_df) <- c("Week","Unvaccinated","Vaccinated")

### dataframe of case by vaccinqation status for plot
week2 <- seq(1,T,1) %>% rep(2)
Type <- c(rep("Unvaccinated",T),rep("Vaccinated",T))
case_df_plot <- cbind(week2,rbind(df_unvac_non,df_vac_non),Type)
colnames(case_df_plot) <- c("Week","Case","Type")