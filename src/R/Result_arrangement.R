result_median <- summary(nuts_fit1)$summary[, "50%"]
result_median <- as.data.frame(result_median)
median_a<-  list()
for(i in 1:500){
  name_a <- paste("a_raw[", i, "]", sep = "")
  median_a$a_raw[i] <- result_median [name_a ,]
  remove(name_a)
}

result_2.5 <- summary(nuts_fit1)$summary[, "2.5%"]
result_2.5 <- as.data.frame(result_2.5)
lower_a<-  list()
for(i in 1:500){
  name_a <- paste("a_raw[", i, "]", sep = "")
  lower_a$a_raw[i] <- result_2.5 [name_a ,]
  remove(name_a)
}

result_97.5 <- summary(nuts_fit1)$summary[, "97.5%"]
result_97.5 <- as.data.frame(result_97.5)
upper_a<-  list()
for(i in 1:500){
  name_a <- paste("a_raw[", i, "]", sep = "")
  upper_a$a_raw[i] <- result_2.5 [name_a ,]
  remove(name_a)
}


result <- cbind(result_median,result_2.5,result_97.5)
colnames(result) <- c("median","lower","upper")
