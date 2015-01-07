install.packages("ggplot2")

library(ggplot2)

# Processing of Ule_EDtable1 sheet 2
 
# Column names within script differ from those within Ule_EDtable1 due to stylistic changes to table during manuscript preparation 

RNASeq_vs_FUS = read.csv("~/RNASeq_vs_FUS.csv", sep=",", head=TRUE)
RNASeq_vs_FUS$ratio.after.before.R <- (as.numeric(RNASeq_vs_FUS$ratio.after.before.R))
RNASeq_vs_FUS$ratio.after.before.F <- (as.numeric(RNASeq_vs_FUS$ratio.after.before.F))
RNASeq_vs_FUS <- RNASeq_vs_FUS[!is.na(RNASeq_vs_FUS$ratio.after.before.R), ]
RNASeq_vs_FUS <- RNASeq_vs_FUS[!is.na(RNASeq_vs_FUS$ratio.after.before.F), ]
RNASeq_vs_FUS <- RNASeq_vs_FUS[!is.infinite(RNASeq_vs_FUS), ]
RNASeq_vs_FUS <- RNASeq_vs_FUS[!is.infinite(RNASeq_vs_FUS$ratio.after.before.R), ]

#Select out non-sig
RNASeq_vs_FUS_non_sig <- subset(RNASeq_vs_FUS, sig == "FALSE")

#Select out sig unproductive going the right way
RNASeq_vs_FUS_sig <- subset(RNASeq_vs_FUS, ratio.after.before.R >"1")
RNASeq_vs_FUS_sig <- subset(RNASeq_vs_FUS_sig, ratio.after.before.F >"1")
RNASeq_vs_FUS_sig <- subset(RNASeq_vs_FUS_sig, sig == "SIG")

#select out sig unproductive going the wrong way and bind to non-sig#
RNASeq_vs_FUS_wrong_way <- subset(RNASeq_vs_FUS_2, ratio.after.before.R <="0")
RNASeq_vs_FUS_wrong_way <- subset(RNASeq_vs_FUS_wrong_way, ratio.after.before.F <="0")
RNASeq_vs_FUS_wrong_way <- subset(RNASeq_vs_FUS_wrong_way, sig == "SIG")
RNASeq_vs_FUS_non_sig_wrong_way <- rbind(RNASeq_vs_FUS_non_sig, RNASeq_vs_FUS_wrong_way)

#select out recursive sites
RNASeq_vs_FUS_RSS_1 <- subset(RNASeq_vs_FUS_sig, Recursive == "YES", )
RNASeq_vs_FUS_RSS_2 <- subset(RNASeq_vs_FUS_sig, Recursive.P == "YES", )
RNASeq_vs_FUS_RSS_3 <- subset(RNASeq_vs_FUS_sig, Recursive.A == "YES", )
RNASeq_vs_FUS_RSS <- rbind(RNASeq_vs_FUS_RSS_1, RNASeq_vs_FUS_RSS_2, RNASeq_vs_FUS_RSS_3)

RNASeq_vs_FUS_non_sig_wrong_way$Extra <- "Non Sig"
RNASeq_vs_FUS_sig$Extra <- "Motif A"
RNASeq_vs_FUS_RSS$Extra <- "Motif B"

RNASeq_vs_FUS_all <- rbind(RNASeq_vs_FUS_non_sig_wrong_way, RNASeq_vs_FUS_sig, RNASeq_vs_FUS_RSS)

ggplot(RNASeq_vs_FUS_all) + geom_point(aes(x = log10(ratio.after.before.R), y = log10(ratio.after.before.F), color=factor(Extra), size=2)) + scale_colour_manual(values=c("black", "red", "grey")) + geom_hline(y=log10(0.894654152), colour="steelblue", linetype = "longdash") + geom_hline(y=log10(1.284250294), colour="steelblue", linetype = "longdash") + geom_vline(x=log10(0.885466603), colour="steelblue", linetype = "longdash") + geom_vline(x=log10(1.198153005), colour="steelblue", linetype = "longdash") + ylim(c(-1.5, 1.5)) + xlim(c(-1.5, 1.5))

