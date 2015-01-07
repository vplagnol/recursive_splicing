install.packages("ggplot2")
install.packages("plyr")

library(ggplot2)
library(plyr)

# Restrict analysis only to those genes were canonical first exon-exon junction detected

Cryptic = read.csv("~/Cryptic_first.csv", sep=",", head=TRUE)

Keep <- Cryptic[grep("Exon-Exon", Cryptic$Junction.type), ]
Keep <- subset(Keep, Counts >= "1")
Cryptic <- Cryptic[which(Cryptic$Gene.ID %in% Keep$Gene.ID),]

# Determine average number of cryptic elements per intron

Cryptic_summary <- ddply(Cryptic, c("Gene.ID", "X1st.Intron.Bin"), summarise, Novel=sum(Junction.type=="Novel"))
Cryptic_summary$X1st.Intron.Bin <- factor(Cryptic_summary$X1st.Intron.Bin, levels = c("0.1-2kb", "2-5kb", "5-10kb", "10-20kb", "20-50kb", "50-100kb", "100kb+"))

ggplot(Cryptic_summary) + geom_boxplot(aes(x=X1st.Intron.Bin, y=Novel), outlier.size = 0, fill="#0072B2") + coord_cartesian(ylim = c(-0.2, 10)) +  scale_y_continuous(breaks=seq(0, 10, 2)) + stat_summary(aes(x=X1st.Intron.Bin, y=Novel), fun.y = "mean", geom = "point",  pch= 18, size= 7, color= "red") 

# Mann Whitney U test with two tails / Wilcoxon Rank Sum Test

G1 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="0.1-2kb")
G2 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="2-5kb")
G3 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="5-10kb")
G4 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="10-20kb")
G5 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="20-50kb")
G6 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="50-100kb")
G7 <- subset(Cryptic_summary, Cryptic_summary$X1st.Intron.Bin=="100kb+")

wilcox.test(G7$Novel, G1$Novel)
wilcox.test(G7$Novel, G2$Novel)
wilcox.test(G7$Novel, G3$Novel)
wilcox.test(G7$Novel, G4$Novel)
wilcox.test(G7$Novel, G5$Novel)
wilcox.test(G7$Novel, G6$Novel)
wilcox.test(G7$Novel, G7$Novel)




