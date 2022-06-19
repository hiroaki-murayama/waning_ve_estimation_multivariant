TM <- 73
X <- seq(1,TM,1)
Y_o <- rep(0,TM)
Y_a <- rep(0,TM)
Y_d <- rep(0,TM)
Y_om <- rep(0,TM)
for(i in 1:TM){
    Y_o[i] <- ve_others(i)
    Y_a[i] <- ve_alpha(i)
    Y_d[i] <- ve_delta(i)
    Y_om[i] <- ve_omicron(i)
}
eps1 <- c(rep(0,l+1),eps)
week <- seq(1,TM,by=1)
week1 <- seq(1,T-1,by=1)
weeks <- seq(1,T+l-1,by=1)
t1 <- TM; t2 <- TM*2; t3 <- TM*3; t4 <- num_data*4 
result[1:TM,] %>% cbind(week,Y_o) -> others_ve
result[(t1+1):(t1+TM),] %>% cbind(week,Y_a) -> alpha_ve
result[(t2+1):(t2+TM),] %>% cbind(week,Y_d) -> delta_ve
result[(t3+1):(t3+TM),] %>% cbind(week,Y_om) -> omicron_ve
result[(t4+1):(t4+T+l-1),] %>% cbind(ncov_vax$week[1:(nrow(A)-1)],ncov_vax$cumulative[1:(nrow(A)-1)],eps1[1:(nrow(A)-1)]) -> overall_ve
colnames(overall_ve) <- c("median","lower","upper","week","vaccine","true_value")
#result[(t3+l+T):(t3+T+l+T-2),] %>% cbind(week1) -> Rit_ve
#result[(t3+2*T+l+32-1):(t3+2*T+T+l+32+-3),] %>% cbind(week1) -> Rjt_ve
#result[(nrow(result)-2*(T-1)):(nrow(result)-T),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> u_pred
#result[(nrow(result)-T+1):(nrow(result)-1),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> v_pred


others_plot <- others_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#2b8cbe", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_o*100),linetype = "dashed", color="#de2d26", size=2, alpha=1)+
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against variant A") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) + 
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100))  

alpha_plot <- alpha_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#2b8cbe", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_a*100),linetype = "dashed", color="#de2d26", size=2, alpha=1)+
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against variant B") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) 

delta_plot <- delta_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#2b8cbe", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_d*100),linetype = "dashed", color="#de2d26", size=2, alpha=1)+
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against variant C") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black")) + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) 

omicron_plot <- omicron_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#2b8cbe", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_om*100),linetype = "dashed", color="#de2d26", size=2, alpha=1)+
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against variant D") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black")) + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) 


scaling_parameter=max(overall_ve$vaccine)/max(overall_ve$upper[!is.na(overall_ve$upper)])
adj = 0.008
options(scipen=1000000)

overall_plot <- overall_ve %>% ggplot() + 
  geom_bar(aes(x=week, y=vaccine/N*100),stat='identity', fill="darkseagreen", width=3) +
  geom_line(aes(x = week, y = median*100), color = "#756bb1", size = 1, alpha=1) +
  geom_ribbon(aes(ymax=lower*100, ymin=upper*100, x=week), fill="#756bb1", alpha = 0.4) +
  labs(x="Week after two doses of vaccination",y="Cumulative number of vaccination with two doses") +
  ggtitle("Overall effect") + theme_bw() +
  theme(text = element_text(size=15, family="sans",color="black"),axis.text.x = element_text(size=15, angle = 90, hjust = 1),axis.text.y = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_date(date_labels="%y/%m/%d",date_breaks  ="7 day", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit=c(0,100),
                     breaks=c(0,25,50,75,100), name="Overall effect of vaccine and\
proportion of vaccinated individuals [%]")


voc_plot <- voc %>% ggplot() +
  geom_bar(aes(x=week,y=Proportion, fill=Variant_type), stat = "identity") +
  scale_x_date(date_labels = "%y/%m/%d", date_breaks = "7 day", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Week", y="Proportion") +
  ggtitle("Proportion of variants") +
  theme_bw() + theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1),axis.text.y = element_text(size=15)) + scale_fill_brewer(palette = "Paired") + 
  theme(text = element_text(size=20, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

case_simulation <- df %>% ggplot() +
  geom_bar(aes(x=Week,y=Case, fill=Type), stat = "identity") +
  scale_x_date(date_labels = "%y/%m/%d", date_breaks = "7 day", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x="Week", y="Incidence") +
  ggtitle("Simulated epidemic curve with vaccination status") + scale_fill_brewer(palette = "Paired") +
  theme_bw() + theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1),axis.text.y = element_text(size=15)) +
  theme(text = element_text(size=20, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
