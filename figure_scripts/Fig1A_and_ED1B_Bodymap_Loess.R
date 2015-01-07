install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")

library(ggplot2)
library(reshape2)
library(plyr)

# Read in DESeq output with which coding genes have been indicated and max gene lengths added

Adipose = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Adipose.csv", sep=",", head=TRUE)
Adrenal = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Adrenal.csv", sep=",", head=TRUE)
Breast = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Breast.csv", sep=",", head=TRUE)
Colon = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Colon.csv", sep=",", head=TRUE)
Heart = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Heart.csv", sep=",", head=TRUE)
Kidney = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Kidney.csv", sep=",", head=TRUE)
Liver = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Liver.csv", sep=",", head=TRUE)
Lung = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Lung.csv", sep=",", head=TRUE)
Lymph = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Lymph_node.csv", sep=",", head=TRUE)
Ovary = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Ovary.csv", sep=",", head=TRUE)
Prostate = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Prostate.csv", sep=",", head=TRUE)
Skeletal = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Skel_muscle.csv", sep=",", head=TRUE)
Testes = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Testes.csv", sep=",", head=TRUE)
Thyroid = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-Thyroid.csv", sep=",", head=TRUE)
White = read.csv("~/tissue-specific-DESeq-nonBinominal-Brain-White_blood_cells.csv", sep=",", head=TRUE)

# Subset relevant data

Adipose <- Adipose[grep("YES", Adipose$coding), ]
Adrenal <- Adrenal[grep("YES", Adrenal$coding), ]
Breast <- Breast[grep("YES", Breast$coding), ]
Colon <- Colon[grep("YES", Colon$coding), ]
Heart <- Heart[grep("YES", Heart$coding), ]
Kidney <- Kidney[grep("YES", Kidney$coding), ]
Liver <- Liver[grep("YES", Liver$coding), ]
Lung <- Lung[grep("YES", Lung$coding), ]
Lymph <- Lymph[grep("YES", Lymph$coding), ]
Ovary <- Ovary[grep("YES", Ovary$coding), ]
Prostate <- Prostate[grep("YES", Prostate$coding), ]
Skeletal <- Skeletal[grep("YES", Skeletal$coding), ]
Testes <- Testes[grep("YES", Testes$coding), ]
Thyroid <- Thyroid[grep("YES", Thyroid$coding), ]
White <- White[grep("YES", White$coding), ]
Brain <- Brain[grep("YES", Brain$coding), ]

Loess <- matrix(nrow=nrow(Adipose), ncol=17)

Loess[, 1] <- Adipose$Length
Loess[, 2] <- log2(Adipose$foldChange)
Loess[, 3] <- log2(Adrenal$foldChange)
Loess[, 4] <- log2(Breast$foldChange)
Loess[, 5] <- log2(Colon$foldChange)
Loess[, 6] <- log2(Heart$foldChange)
Loess[, 7] <- log2(Kidney$foldChange)
Loess[, 8] <- log2(Liver$foldChange)
Loess[, 9] <- log2(Lung$foldChange)
Loess[, 10] <- log2(Lymph$foldChange)
Loess[, 11] <- log2(Ovary$foldChange)
Loess[, 12] <- log2(Prostate$foldChange)
Loess[, 13] <- log2(Skeletal$foldChange)
Loess[, 14] <- log2(Testes$foldChange)
Loess[, 15] <- log2(Thyroid$foldChange)
Loess[, 16] <- log2(White$foldChange)
Loess[, 17] <- log2(1)

colnames(Loess) <- c("Length", "Adipose", "Adrenal", "Breast", "Colon", "Heart", "Kidney", "Liver", "Lung", "Lymph", "Ovary", "Prostate", "Skeletal", "Testes", "Thyroid", "White", "Brain") 

Loess <- as.data.frame(Loess)

Loess <- melt(Loess, id.vars="Length")
Loess$Length <- (as.numeric(Loess$Length))
Loess$value <- (as.numeric(Loess$value))
Loess <- Loess[!is.na(Loess$value), ]
Loess <- Loess[!is.infinite(Loess$value), ]
Loess <- Loess[!is.na(Loess$Length), ]
Loess <- Loess[!is.infinite(Loess$Length), ]

ggplot(Loess, aes(x=log10(Length), y=value, colour=factor(variable))) + geom_smooth(method = "loess", aes(x=log10(Length), y=value), se = FALSE, size=1.2) + ylab("Log2(Fold Expression vs. Brain)") + xlab("Gene length") + geom_vline(x=log10(150000), colour="steelblue", linetype = "longdash") 

