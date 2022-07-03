### check fitness by comparing observed case data with simulated cases from posterior samples
unvaccinated = (ncovdata$unvax_case[(2+l):(T+l)] - ncovdata$exclusion[(2+l):(T+l)]) + ncovdata$unknown[(2+l):(T+l)] * ((ncovdata$unvax_case[(2+l):(T+l)] - ncovdata$exclusion[(2+l):(T+l)])/((ncovdata$unvax_case[(2+l):(T+l)] - ncovdata$exclusion[(2+l):(T+l)])+ncovdata$vax_case[(2+l):(T+l)])) 
vaccinated = ncovdata$vax_case[(2+l):(T+l)] + ncovdata$unknown[(2+l):(T+l)] * ((ncovdata$vax_case[(2+l):(T+l)] - ncovdata$exclusion[(2+l):(T+l)])/((ncovdata$unvax_case[(2+l):(T+l)] - ncovdata$exclusion[(2+l):(T+l)])+ncovdata$vax_case[(2+l):(T+l)]))
u_pred %<>% cbind(unvaccinated); colnames(u_pred) <- c("Median","Lower","Upper","Week","Observed"); 
v_pred %<>% cbind(vaccinated); colnames(v_pred) <- c("Median","Lower","Upper","Week","Observed")

week <- as.data.frame(c(u_pred$Week,u_pred$Week))
type <- as.data.frame(c(rep("Observed",T-1),rep("Predicted",T-1)))
case <- as.data.frame(c(u_pred$Observed,u_pred$Median))
lower <- as.data.frame(c(rep(NaN,T-1),u_pred$Lower))
upper <- as.data.frame(c(rep(NaN,T-1),u_pred$Upper))
u_pred <- cbind(week,type,case,lower,upper)
colnames(u_pred) <- c("Week","Type","case","lower","upper")
u_pred <- as.data.frame(u_pred)

week <- as.data.frame(c(v_pred$Week,v_pred$Week))
type <- as.data.frame(c(rep("Observed",T-1),rep("Predicted",T-1)))
case <- as.data.frame(c(v_pred$Observed,v_pred$Median))
lower <- as.data.frame(c(rep(NaN,T-1),v_pred$Lower))
upper <- as.data.frame(c(rep(NaN,T-1),v_pred$Upper))
v_pred <- cbind(week,type,case,lower,upper)
colnames(v_pred) <- c("Week","Type","case","lower","upper")
u_pred <- as.data.frame(u_pred)

u_param_plot <- ggplot(u_pred, aes(x=Week-2, y=case, fill=Type)) + scale_fill_brewer(palette = "Paired") + labs(x="Date", y="Unvaccinated incidence") + ggtitle("Comparison between observed incidence and predicted estimates") +
  geom_bar(stat='identity', position='dodge') +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymin=lower,ymax=upper),width=1,position=position_dodge(4.5)) +
  theme(text = element_text(size=20, family="sans",color="black"),
        axis.text = element_text(size=18, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") + #geom_text(aes(label = age), size = 2, hjust = 0.5, vjust = 3, position = "stack")  + 
  scale_x_date(date_labels="%d/%m/%y",date_breaks  ="7 day", expand = c(0, 0)) + guides(fill=guide_legend(""))

v_param_plot <- ggplot(v_pred, aes(x=Week-2, y=case, fill=Type)) + scale_fill_brewer(palette = "Paired") + labs(x="Date", y="Vaccinated incidence") + ggtitle("Comparison between observed incidence and predicted estimates") +
  geom_bar(stat='identity', position='dodge') +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymin=lower,ymax=upper),width=1,position=position_dodge(4.5)) +
  theme(text = element_text(size=20, family="sans",color="black"),
        axis.text = element_text(size=18, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") + #geom_text(aes(label = age), size = 2, hjust = 0.5, vjust = 3, position = "stack")  + 
  scale_x_date(date_labels="%d/%m/%y",date_breaks  ="7 day", expand = c(0, 0)) + guides(fill=guide_legend(""))