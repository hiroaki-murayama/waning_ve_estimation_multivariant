TM <- 73
X <- seq(1,TM,1)
Y_o <- rep(0,TM)
Y_a <- rep(0,TM)
Y_d <- rep(0,TM)
Y_om <- rep(0,TM)
for(i in 1:TM){
    Y_o[i] <- ve_A(i)
    Y_a[i] <- ve_B(i)
    Y_d[i] <- ve_C(i)
    Y_om[i] <- ve_D(i)
}
week <- seq(1,TM,by=1)
week1 <- seq(1,(T+l-1),by=1)
weeks <- seq(1,T+l-1,by=1)
t1 <- TM; t2 <- TM*2; t3 <- TM*3; t4 <- num_data*4 
result[1:TM,] %>% cbind(week,Y_o) -> others_ve
result[(t1+1):(t1+TM),] %>% cbind(week,Y_a) -> alpha_ve
result[(t2+1):(t2+TM),] %>% cbind(week,Y_d) -> delta_ve
result[(t3+1):(t3+TM),] %>% cbind(week,Y_om) -> omicron_ve
result[(t4+1):(t4+T+l-1),] %>% cbind(week1,ncov_vax$cumulative[1:(T+l-1)]) -> overall_ve
colnames(overall_ve) <- c("median","lower","upper","week","vaccine")
N = 15000000

others_plot <- others_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#0092A7", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#0092A7", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_o*100),linetype = "dashed", color="#6F8793", size=2, alpha=1)+
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
  ggtitle("Waning vaccine effectiveness against variant A") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) + 
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100))  

alpha_plot <- alpha_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#0092A7", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#0092A7", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_a*100),linetype = "dashed", color="#6F8793", size=2, alpha=1)+
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
  ggtitle("Waning vaccine effectiveness against variant B") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, TM, 4), limits=c(1, TM)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) 

delta_plot <- delta_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#0092A7", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#0092A7", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_d*100),linetype = "dashed", color="#6F8793", size=2, alpha=1)+
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
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
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#0092A7", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), color="#0092A7", size = 1, alpha=1) +
  geom_line(aes(x = week, y = Y_om*100),linetype = "dashed", color="#6F8793", size=2, alpha=1)+
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
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



overall_plot <- overall_ve %>% ggplot() + 
  geom_bar(aes(x=week, y=vaccine/N*100),stat='identity', fill="darkseagreen", width=0.5) +
  geom_line(aes(x = week, y = median*100), color = "#9D7F4E", size = 1, alpha=1) +
  geom_ribbon(aes(ymax=lower*100, ymin=upper*100, x=week), fill="#9D7F4E", alpha = 0.4) +
  labs(x="Week since completion of required dose(s)",y="Cumulative number of vaccination with two doses") +
  theme_bw() +
  theme(text = element_text(size=15, family="sans",color="black"),axis.text.x = element_text(size=15, hjust = 1),axis.text.y = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(1, T+l-1, 1), limits=c(0, T+l), expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limit=c(0,100),
                     breaks=c(0,25,50,75,100), name="Population-level protection\
Proportion vaccinated [%]")


voc_plot <- case_voc %>% ggplot() +
  geom_bar(aes(x=Week,y=Proportion, fill=Variant), stat = "identity") +
  scale_x_continuous(breaks = seq(1, 73, 4), limits=c(0, 74), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x="Week", y="Proportion") +
  ggtitle("Assumed proportion of variants") +
  theme_bw() + theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15)) + scale_fill_brewer(palette = "Paired") + 
    theme(text = element_text(size=20, family="sans",color="black"),
          axis.text = element_text(size=10, family="sans",color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(fill=guide_legend(""))

case_simulation <- case_df_plot %>% ggplot() +
  geom_bar(aes(x=Week,y=Case, fill=Type), stat = "identity", width=0.5) +
  scale_x_continuous(breaks = seq(1, 73, 4), limits=c(0, 74), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x="Week", y="Incidence") +
  ggtitle("Simulated epidemic curve with vaccination status") + scale_fill_brewer(palette = "Paired") +
  theme_bw() + theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15)) +
    theme(text = element_text(size=20, family="sans",color="black"),
          axis.text = element_text(size=10, family="sans",color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(fill=guide_legend(""))
