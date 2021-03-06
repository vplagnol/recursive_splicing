library(ggplot2)
library(reshape2)
library(plyr)


Adipose = read.csv("data/RNASeq_bodymap/Adipose_vs_brain.csv.gz", sep=",", head=TRUE)
Adrenal = read.csv("data/RNASeq_bodymap/Adrenal_vs_brain.csv.gz", sep=",", head=TRUE)
Breast = read.csv("data/RNASeq_bodymap/Breast_vs_brain.csv.gz", sep=",", head=TRUE)
Colon = read.csv("data/RNASeq_bodymap/Colon_vs_brain.csv.gz", sep=",", head=TRUE)
Heart = read.csv("data/RNASeq_bodymap/Heart_vs_brain.csv.gz", sep=",", head=TRUE)
Kidney = read.csv("data/RNASeq_bodymap/Kidney_vs_brain.csv.gz", sep=",", head=TRUE)
Liver = read.csv("data/RNASeq_bodymap/Liver_vs_brain.csv.gz", sep=",", head=TRUE)
Lung = read.csv("data/RNASeq_bodymap/Lung_vs_brain.csv.gz", sep=",", head=TRUE)
Lymph = read.csv("data/RNASeq_bodymap/Lymph_vs_brain.csv.gz", sep=",", head=TRUE)
Ovary = read.csv("data/RNASeq_bodymap/Ovary_vs_brain.csv.gz", sep=",", head=TRUE)
Prostate = read.csv("data/RNASeq_bodymap/Prostate_vs_brain.csv.gz", sep=",", head=TRUE)
Skeletal = read.csv("data/RNASeq_bodymap/Skel_muscle_vs_brain.csv.gz", sep=",", head=TRUE)
Testes = read.csv("data/RNASeq_bodymap/Testes_vs_brain.csv.gz", sep=",", head=TRUE)
Thyroid = read.csv("data/RNASeq_bodymap/Thyroid_vs_brain.csv.gz", sep=",", head=TRUE)
White = read.csv("data/RNASeq_bodymap/White_vs_brain.csv.gz", sep=",", head=TRUE)

message('Done parsing the input files')


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

p <- ggplot(Loess, aes(x=log10(Length), y=value, colour=factor(variable))) + geom_smooth(method = "loess", aes(x=log10(Length), y=value), se = FALSE, size=1.2) + ylab("Log2(Fold Expression vs. Brain)") + xlab("Gene length") + geom_vline(x=log10(150000), colour="steelblue", linetype = "longdash") 


ggsave(plot = p, filename = 'fig/bodymap_loess.pdf')