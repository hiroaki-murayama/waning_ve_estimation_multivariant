#####outcome arrangement#####
week <- seq(1,num_data,by=1)
week1 <- seq(1,T-1,by=1)
t1 <- (num_data); t2 <- (num_data)*2; t3 <- (num_data)*3
result[(1):(num_data),] %>% cbind(week,lower_ref,median_ref,upper_ref,lower_ref1,median_ref1,upper_ref1) -> delta_ve
result[(t1+1):(t1+num_data),] %>% cbind(week,lower_ref_om,median_ref_om,upper_ref_om) -> omicron_ve
result[(t2+1):(t2+num_data-1),] %>% cbind(ncov_vax$week[1:(num_data-1)],ncov_vax$cumulative[1:(num_data-1)]) -> overall_ve
colnames(overall_ve) <- c("median","lower","upper","week","vaccine")
result[(nrow(result)-2*(T-1)):(nrow(result)-T),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> u_pred
result[(nrow(result)-T+1):(nrow(result)-1),] %>% cbind(ncovdata$week[(l+1):(T+l-1)]) -> v_pred

#####plot#####
delta_plot <- delta_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower_ref, ymax = upper_ref), fill = "#f03b20", alpha = 0.4) +
  geom_line(aes(x = week, y = median_ref, colour="British Columbia"), size = 1, alpha=1) +
  geom_ribbon(aes(x = week, ymin = lower_ref1, ymax = upper_ref1), fill = "#756bb1", alpha = 0.4) +
  geom_line(aes(x = week, y = median_ref1, colour="Quebec"), size = 1, alpha=1) +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), colour="#2b8cbe", size = 1, alpha=1) +
  scale_colour_manual(name=c(""), values = c("Estimated waning vaccine effectiveness"="#2b8cbe", "British Columbia"="#f03b20", "Quebec"="#756bb1")) +
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
  ggtitle("Waning vaccine effectiveness against Delta in Japan") +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, num_data, 4), limits=c(1, num_data)) + 
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) + theme(legend.position=c(.3,.2)) + guides(fill=guide_legend(""))

omicron_plot <- omicron_ve %>% ggplot() +
  geom_ribbon(aes(x = week, ymin = lower_ref_om, ymax = upper_ref_om), fill = "#31a354", alpha = 0.4) +
  geom_line(aes(x = week, y = median_ref_om, colour="Japan (NIID)"), size = 1, alpha=1) +
  geom_ribbon(aes(x = week, ymin = lower*100, ymax = upper*100), fill = "#2b8cbe", alpha = 0.4) +
  geom_line(aes(x = week, y = median*100), colour="#2b8cbe", size = 1, alpha=1) +
  labs(x="Week since completion of required dose(s)",y="Vaccine effectiveness [%]") +
  ggtitle("Waning vaccine effectiveness against Omicron in Japan") +
  scale_colour_manual(name=c(""), values = c("Estimated waning vaccine effectiveness"="#2b8cbe", "Japan (NIID)"="#31a354")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1),axis.text.y = element_text(size=15),
        text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=10, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right") + scale_x_continuous(breaks = seq(1, num_data, 4), limits=c(1, num_data)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits=c(0, 100)) + theme(legend.position=c(.3,.2)) + guides(fill=guide_legend("")) 

overall_plot <- overall_ve %>% ggplot() + 
  geom_bar(aes(x=week, y=vaccine*100),stat='identity', fill="darkseagreen", width=3) +
  geom_line(aes(x = week, y = median*100), color = "#9D7F4E", size = 1, alpha=1) +
  geom_ribbon(aes(ymax=lower*100, ymin=upper*100, x=week), fill="#9D7F4E", alpha = 0.4) +
  labs(x="Week since completion of required dose(s)",y="Cumulative number of vaccination with required dose(s)") +
  theme_bw() +
  theme(text = element_text(size=15, family="sans",color="black"),axis.text.x = element_text(size=15, angle = 90, hjust = 1),axis.text.y = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_date(date_labels="%d/%m/%y",date_breaks  ="7 day", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limit=c(0,100),
                     breaks=c(0,25,50,75,100), name="Population-level protection\
Proportion vaccinated [%]") + guides(fill=guide_legend(""))
