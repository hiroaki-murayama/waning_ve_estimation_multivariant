#####outcome arrangement#####
week <- seq(1,30,by=1)
week1 <- seq(1,T-1,by=1)
weeks <- seq(1,T+l+32-1,by=1)
t1 <- (num_data); t2 <- (num_data)*2; t3 <- (num_data)*3
result[1:30,] %>% cbind(week,lower_ref,median_ref,upper_ref) -> others_ve
result[(t1+1):(t1+30),] %>% cbind(week,lower_ref,median_ref,upper_ref) -> alpha_ve
result[(t2+1):(t2+30),] %>% cbind(week,lower_ref,median_ref,upper_ref) -> delta_ve
result[(t3+1):(t3+num_data-1),] %>% cbind(ncov_vax$week[1:(num_data-1)],ncov_vax$cumulative[1:(num_data-1)]) -> overall_ve
colnames(overall_ve) <- c("median","lower","upper","week","vaccine")
result[(nrow(result)-2*(T-1)):(nrow(result)-T),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> u_pred
result[(nrow(result)-T+1):(nrow(result)-1),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> v_pred

#####plot#####
others_plot <- others_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), colour="#2b8cbe", size = 1, alpha=1) +
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against Others in Ontario") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, 30, 4), limits=c(1, 30)) + 
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100))  

alpha_plot <- alpha_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), colour="#2b8cbe", size = 1, alpha=1) +
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against Alpha in Ontario") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, 30, 4), limits=c(1, 30)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) 

delta_plot <- delta_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), colour="#2b8cbe", size = 1, alpha=1) +
  geom_ribbon(aes(x = week, ymin = lower_ref, ymax = upper_ref), fill = "#f03b20", alpha = 0.4) +
  geom_line(aes(x = week, y = median_ref, colour="British Columbia"), size = 1, alpha=1) +
  geom_ribbon(aes(x = week, ymin = lower_ref1, ymax = upper_ref1), fill = "#756bb1", alpha = 0.4) +
  geom_line(aes(x = week, y = median_ref1, colour="Quebec"), size = 1, alpha=1) +
  labs(x="Week after two doses of vaccination",y="Direct effect of vaccine [%]") +
  ggtitle("Waning vaccine effectiveness against Delta in Ontario") +
  scale_colour_manual(name="Reference", values = c("British Columbia"="#f03b20", "Quebec"="#756bb1")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1), axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",colour="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = seq(1, 30, 4), limits=c(1, 30)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) +theme(legend.position=c(.2,.3),legend.background = element_rect(fill = "white", colour = "black")) 


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
