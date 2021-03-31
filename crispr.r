setwd("C:/Users/sirui.zhou/work/Laetitia.crispr")
library(data.table)
library(dplyr)
count <- fread("new.trim40-51.count.txt", header=T)
nrow(count)
hist(count$plasmidnew)

#count %>% 
#top10 <- arrange(desc(plasmidnew)) %>% filter(plasmidnew > quantile(plasmidnew, .9))
#filter(plasmidnew > quantile(plasmidnew, .1))

count2 <- count[ which(count$plasmidnew != 0, ), ]

count3 <- count[ which(count2$plasmidnew < 20000 & count2$plasmidnew > 0), ]
count4 <- count[ which(count2$plasmidnew < 10000 & count2$plasmidnew > 0), ]
count5 <- count[ which(count2$plasmidnew < 10000 & count2$plasmidnew > 2), ]
nrow(count2)
nrow(count3)
nrow(count4)
nrow(count5)

###calculate skew ratio###

q10 <- quantile(count$plasmidnew, .1)
q90 <- quantile(count$plasmidnew, .9)

mean(q90)
mean(q10)

top10 <- count5 %>% top_n(1005,plasmidnew)
bot10 <- count5 %>% top_n(-1005,plasmidnew)


nrow(top10)
nrow(bot10)


top10$quant <- "top10"
bot10$quant <- "bot10"



dat <- rbind(top10, bot10)

boxplot(top10$plasmidnew, bot10$plasmidnew, main="top10 vs bottom10 after QC",
          names=c("top10", "bottom10"), col=c("red","blue") )
mean(top10$plasmidnew) - mean(bot10$plasmidnew)


library(ggplot2)

count2$log2 <- log2(count2$plasmidnew)

q10 <- quantile(count5$plasmidnew, .1)
q90 <- quantile(count5$plasmidnew, .9)

mean(q90)
mean(q10)
mean(q90)/mean(q10)

p <- ggplot(count2, aes(x=log2)) + 
  geom_density() +
  labs(title="Counts Density", 
       #subtitle = "No zero counts, remove sgRNA counts > 10000 AND < 2, sgRNA#=10047",  
       x = "Log2 Counts" )
p+scale_color_grey() + theme_classic() + geom_vline(aes(xintercept=mean(log2)),
                                                    color="blue", linetype="dashed", size=1)

p <- ggplot(dat, aes(x=plasmidnew, color=quant)) + 
  geom_density() +
  labs(title="Counts Density", subtitle = "No zero counts, remove sgRNA counts > 10000 AND < 2, sgRNA#=10047",  x = "Counts" )
library(plyr)
mu <- ddply(dat, "quant", summarise, grp.mean=mean(plasmidnew))
p + theme_classic() + geom_vline(data=mu, aes(xintercept=grp.mean, color=quant),
           linetype="dashed")
