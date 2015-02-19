library(ggplot2)
library(reshape2)
library(plyr)
library(plotrix)
library(gridExtra)


# Read in DESeq output with which coding genes have been indicated and max gene lengths added
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

# Process data for Adipose

Adipose <- Adipose[grep("YES", Adipose$coding), ]
Adipose <- subset(Adipose, baseMeanA>0, baseMeanB>0)
Adipose <- subset(Adipose, foldChange>0)
Adipose2 <- melt(Adipose, id.vars="Length", "foldChange")
Adipose2$value <- log2(as.numeric(Adipose2$value))
Adipose2$Length <- (as.numeric(Adipose2$Length))
Adipose2 <- Adipose2[!is.na(Adipose2$Length), ]
Adipose2 <- Adipose2[!is.na(Adipose2$value), ]
Adipose2 <- Adipose2[!is.infinite(Adipose2$Length), ]
Adipose2 <- Adipose2[!is.infinite(Adipose2$value), ]
Adipose2$Type <- "Other"

# Subset recursive genes

Adipose3 <- subset(Adipose, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Adipose3 <- melt(Adipose3, id.vars="Length", "foldChange")
Adipose3$value <- log2(as.numeric(Adipose3$value))
Adipose3$Length <- (as.numeric(Adipose3$Length))
Adipose3 <- Adipose3[!is.na(Adipose3$Length), ]
Adipose3 <- Adipose3[!is.na(Adipose3$value), ]
Adipose3 <- Adipose3[!is.infinite(Adipose3$Length), ]
Adipose3 <- Adipose3[!is.infinite(Adipose3$value), ]
Adipose3$Type <- "RSS"

# Subset Titin
Adipose4 <- subset(Adipose, id=="ENSG00000155657")
Adipose4 <- melt(Adipose4, id.vars="Length", "foldChange")
Adipose4$value <- log2(as.numeric(Adipose4$value))
Adipose4$Length <- (as.numeric(Adipose4$Length))
Adipose4 <- Adipose4[!is.na(Adipose4$Length), ]
Adipose4 <- Adipose4[!is.na(Adipose4$value), ]
Adipose4 <- Adipose4[!is.infinite(Adipose4$Length), ]
Adipose4 <- Adipose4[!is.infinite(Adipose4$value), ]
Adipose4$Type <- "TTN"

# Subset Dystrophin
Adipose5 <- subset(Adipose, id=="ENSG00000198947")
Adipose5 <- melt(Adipose5, id.vars="Length", "foldChange")
Adipose5$value <- log2(as.numeric(Adipose5$value))
Adipose5$Length <- (as.numeric(Adipose5$Length))
Adipose5 <- Adipose5[!is.na(Adipose5$Length), ]
Adipose5 <- Adipose5[!is.na(Adipose5$value), ]
Adipose5 <- Adipose5[!is.infinite(Adipose5$Length), ]
Adipose5 <- Adipose5[!is.infinite(Adipose5$value), ]
Adipose5$Type <- "DMD"

# Merge and plot

Adipose6 <- rbind(Adipose2, Adipose3, Adipose4, Adipose5) 
Adipose_scatter <- ggplot(Adipose6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Adipose2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Adipose") + theme(legend.position="none")

# Process data for Adrenal

Adrenal <- Adrenal[grep("YES", Adrenal$coding), ]
Adrenal <- subset(Adrenal, baseMeanA>0, baseMeanB>0)
Adrenal <- subset(Adrenal, foldChange>0)
Adrenal2 <- melt(Adrenal, id.vars="Length", "foldChange")
Adrenal2$value <- log2(as.numeric(Adrenal2$value))
Adrenal2$Length <- (as.numeric(Adrenal2$Length))
Adrenal2 <- Adrenal2[!is.na(Adrenal2$Length), ]
Adrenal2 <- Adrenal2[!is.na(Adrenal2$value), ]
Adrenal2 <- Adrenal2[!is.infinite(Adrenal2$Length), ]
Adrenal2 <- Adrenal2[!is.infinite(Adrenal2$value), ]
Adrenal2$Type <- "Other"

# Subset recursive genes

Adrenal3 <- subset(Adrenal, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Adrenal3 <- melt(Adrenal3, id.vars="Length", "foldChange")
Adrenal3$value <- log2(as.numeric(Adrenal3$value))
Adrenal3$Length <- (as.numeric(Adrenal3$Length))
Adrenal3 <- Adrenal3[!is.na(Adrenal3$Length), ]
Adrenal3 <- Adrenal3[!is.na(Adrenal3$value), ]
Adrenal3 <- Adrenal3[!is.infinite(Adrenal3$Length), ]
Adrenal3 <- Adrenal3[!is.infinite(Adrenal3$value), ]
Adrenal3$Type <- "RSS"

# Subset Titin
Adrenal4 <- subset(Adrenal, id=="ENSG00000155657")
Adrenal4 <- melt(Adrenal4, id.vars="Length", "foldChange")
Adrenal4$value <- log2(as.numeric(Adrenal4$value))
Adrenal4$Length <- (as.numeric(Adrenal4$Length))
Adrenal4 <- Adrenal4[!is.na(Adrenal4$Length), ]
Adrenal4 <- Adrenal4[!is.na(Adrenal4$value), ]
Adrenal4 <- Adrenal4[!is.infinite(Adrenal4$Length), ]
Adrenal4 <- Adrenal4[!is.infinite(Adrenal4$value), ]
Adrenal4$Type <- "TTN"

# Subset Dystrophin
Adrenal5 <- subset(Adrenal, id=="ENSG00000198947")
Adrenal5 <- melt(Adrenal5, id.vars="Length", "foldChange")
Adrenal5$value <- log2(as.numeric(Adrenal5$value))
Adrenal5$Length <- (as.numeric(Adrenal5$Length))
Adrenal5 <- Adrenal5[!is.na(Adrenal5$Length), ]
Adrenal5 <- Adrenal5[!is.na(Adrenal5$value), ]
Adrenal5 <- Adrenal5[!is.infinite(Adrenal5$Length), ]
Adrenal5 <- Adrenal5[!is.infinite(Adrenal5$value), ]
Adrenal5$Type <- "DMD"

# Merge and plot

Adrenal6 <- rbind(Adrenal2, Adrenal3, Adrenal4, Adrenal5) 
Adrenal_scatter <- ggplot(Adrenal6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Adrenal2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Adrenal") + theme(legend.position="none")


# Process data for Breast

Breast <- Breast[grep("YES", Breast$coding), ]
Breast <- subset(Breast, baseMeanA>0, baseMeanB>0)
Breast <- subset(Breast, foldChange>0)
Breast2 <- melt(Breast, id.vars="Length", "foldChange")
Breast2$value <- log2(as.numeric(Breast2$value))
Breast2$Length <- (as.numeric(Breast2$Length))
Breast2 <- Breast2[!is.na(Breast2$Length), ]
Breast2 <- Breast2[!is.na(Breast2$value), ]
Breast2 <- Breast2[!is.infinite(Breast2$Length), ]
Breast2 <- Breast2[!is.infinite(Breast2$value), ]
Breast2$Type <- "Other"

# Subset recursive genes

Breast3 <- subset(Breast, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Breast3 <- melt(Breast3, id.vars="Length", "foldChange")
Breast3$value <- log2(as.numeric(Breast3$value))
Breast3$Length <- (as.numeric(Breast3$Length))
Breast3 <- Breast3[!is.na(Breast3$Length), ]
Breast3 <- Breast3[!is.na(Breast3$value), ]
Breast3 <- Breast3[!is.infinite(Breast3$Length), ]
Breast3 <- Breast3[!is.infinite(Breast3$value), ]
Breast3$Type <- "RSS"

# Subset Titin
Breast4 <- subset(Breast, id=="ENSG00000155657")
Breast4 <- melt(Breast4, id.vars="Length", "foldChange")
Breast4$value <- log2(as.numeric(Breast4$value))
Breast4$Length <- (as.numeric(Breast4$Length))
Breast4 <- Breast4[!is.na(Breast4$Length), ]
Breast4 <- Breast4[!is.na(Breast4$value), ]
Breast4 <- Breast4[!is.infinite(Breast4$Length), ]
Breast4 <- Breast4[!is.infinite(Breast4$value), ]
Breast4$Type <- "TTN"

# Subset Dystrophin
Breast5 <- subset(Breast, id=="ENSG00000198947")
Breast5 <- melt(Breast5, id.vars="Length", "foldChange")
Breast5$value <- log2(as.numeric(Breast5$value))
Breast5$Length <- (as.numeric(Breast5$Length))
Breast5 <- Breast5[!is.na(Breast5$Length), ]
Breast5 <- Breast5[!is.na(Breast5$value), ]
Breast5 <- Breast5[!is.infinite(Breast5$Length), ]
Breast5 <- Breast5[!is.infinite(Breast5$value), ]
Breast5$Type <- "DMD"

# Merge and plot

Breast6 <- rbind(Breast2, Breast3, Breast4, Breast5) 
Breast_scatter <- ggplot(Breast6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Breast2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Breast") + theme(legend.position="none")


# Process data for Colon

Colon <- Colon[grep("YES", Colon$coding), ]
Colon <- subset(Colon, baseMeanA>0, baseMeanB>0)
Colon <- subset(Colon, foldChange>0)
Colon2 <- melt(Colon, id.vars="Length", "foldChange")
Colon2$value <- log2(as.numeric(Colon2$value))
Colon2$Length <- (as.numeric(Colon2$Length))
Colon2 <- Colon2[!is.na(Colon2$Length), ]
Colon2 <- Colon2[!is.na(Colon2$value), ]
Colon2 <- Colon2[!is.infinite(Colon2$Length), ]
Colon2 <- Colon2[!is.infinite(Colon2$value), ]
Colon2$Type <- "Other"

# Subset recursive genes

Colon3 <- subset(Colon, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Colon3 <- melt(Colon3, id.vars="Length", "foldChange")
Colon3$value <- log2(as.numeric(Colon3$value))
Colon3$Length <- (as.numeric(Colon3$Length))
Colon3 <- Colon3[!is.na(Colon3$Length), ]
Colon3 <- Colon3[!is.na(Colon3$value), ]
Colon3 <- Colon3[!is.infinite(Colon3$Length), ]
Colon3 <- Colon3[!is.infinite(Colon3$value), ]
Colon3$Type <- "RSS"

# Subset Titin
Colon4 <- subset(Colon, id=="ENSG00000155657")
Colon4 <- melt(Colon4, id.vars="Length", "foldChange")
Colon4$value <- log2(as.numeric(Colon4$value))
Colon4$Length <- (as.numeric(Colon4$Length))
Colon4 <- Colon4[!is.na(Colon4$Length), ]
Colon4 <- Colon4[!is.na(Colon4$value), ]
Colon4 <- Colon4[!is.infinite(Colon4$Length), ]
Colon4 <- Colon4[!is.infinite(Colon4$value), ]
Colon4$Type <- "TTN"

# Subset Dystrophin
Colon5 <- subset(Colon, id=="ENSG00000198947")
Colon5 <- melt(Colon5, id.vars="Length", "foldChange")
Colon5$value <- log2(as.numeric(Colon5$value))
Colon5$Length <- (as.numeric(Colon5$Length))
Colon5 <- Colon5[!is.na(Colon5$Length), ]
Colon5 <- Colon5[!is.na(Colon5$value), ]
Colon5 <- Colon5[!is.infinite(Colon5$Length), ]
Colon5 <- Colon5[!is.infinite(Colon5$value), ]
Colon5$Type <- "DMD"

# Merge and plot

Colon6 <- rbind(Colon2, Colon3, Colon4, Colon5) 
Colon_scatter <- ggplot(Colon6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Colon2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Colon") + theme(legend.position="none")


# Process data for Heart

Heart <- Heart[grep("YES", Heart$coding), ]
Heart <- subset(Heart, baseMeanA>0, baseMeanB>0)
Heart <- subset(Heart, foldChange>0)
Heart2 <- melt(Heart, id.vars="Length", "foldChange")
Heart2$value <- log2(as.numeric(Heart2$value))
Heart2$Length <- (as.numeric(Heart2$Length))
Heart2 <- Heart2[!is.na(Heart2$Length), ]
Heart2 <- Heart2[!is.na(Heart2$value), ]
Heart2 <- Heart2[!is.infinite(Heart2$Length), ]
Heart2 <- Heart2[!is.infinite(Heart2$value), ]
Heart2$Type <- "Other"

# Subset recursive genes

Heart3 <- subset(Heart, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Heart3 <- melt(Heart3, id.vars="Length", "foldChange")
Heart3$value <- log2(as.numeric(Heart3$value))
Heart3$Length <- (as.numeric(Heart3$Length))
Heart3 <- Heart3[!is.na(Heart3$Length), ]
Heart3 <- Heart3[!is.na(Heart3$value), ]
Heart3 <- Heart3[!is.infinite(Heart3$Length), ]
Heart3 <- Heart3[!is.infinite(Heart3$value), ]
Heart3$Type <- "RSS"

# Subset Titin
Heart4 <- subset(Heart, id=="ENSG00000155657")
Heart4 <- melt(Heart4, id.vars="Length", "foldChange")
Heart4$value <- log2(as.numeric(Heart4$value))
Heart4$Length <- (as.numeric(Heart4$Length))
Heart4 <- Heart4[!is.na(Heart4$Length), ]
Heart4 <- Heart4[!is.na(Heart4$value), ]
Heart4 <- Heart4[!is.infinite(Heart4$Length), ]
Heart4 <- Heart4[!is.infinite(Heart4$value), ]
Heart4$Type <- "TTN"

# Subset Dystrophin
Heart5 <- subset(Heart, id=="ENSG00000198947")
Heart5 <- melt(Heart5, id.vars="Length", "foldChange")
Heart5$value <- log2(as.numeric(Heart5$value))
Heart5$Length <- (as.numeric(Heart5$Length))
Heart5 <- Heart5[!is.na(Heart5$Length), ]
Heart5 <- Heart5[!is.na(Heart5$value), ]
Heart5 <- Heart5[!is.infinite(Heart5$Length), ]
Heart5 <- Heart5[!is.infinite(Heart5$value), ]
Heart5$Type <- "DMD"

# Merge and plot

Heart6 <- rbind(Heart2, Heart3, Heart4, Heart5) 
Heart_scatter <- ggplot(Heart6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Heart2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Heart") + theme(legend.position="none")


# Process data for Kidney

Kidney <- Kidney[grep("YES", Kidney$coding), ]
Kidney <- subset(Kidney, baseMeanA>0, baseMeanB>0)
Kidney <- subset(Kidney, foldChange>0)
Kidney2 <- melt(Kidney, id.vars="Length", "foldChange")
Kidney2$value <- log2(as.numeric(Kidney2$value))
Kidney2$Length <- (as.numeric(Kidney2$Length))
Kidney2 <- Kidney2[!is.na(Kidney2$Length), ]
Kidney2 <- Kidney2[!is.na(Kidney2$value), ]
Kidney2 <- Kidney2[!is.infinite(Kidney2$Length), ]
Kidney2 <- Kidney2[!is.infinite(Kidney2$value), ]
Kidney2$Type <- "Other"

# Subset recursive genes

Kidney3 <- subset(Kidney, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Kidney3 <- melt(Kidney3, id.vars="Length", "foldChange")
Kidney3$value <- log2(as.numeric(Kidney3$value))
Kidney3$Length <- (as.numeric(Kidney3$Length))
Kidney3 <- Kidney3[!is.na(Kidney3$Length), ]
Kidney3 <- Kidney3[!is.na(Kidney3$value), ]
Kidney3 <- Kidney3[!is.infinite(Kidney3$Length), ]
Kidney3 <- Kidney3[!is.infinite(Kidney3$value), ]
Kidney3$Type <- "RSS"

# Subset Titin
Kidney4 <- subset(Kidney, id=="ENSG00000155657")
Kidney4 <- melt(Kidney4, id.vars="Length", "foldChange")
Kidney4$value <- log2(as.numeric(Kidney4$value))
Kidney4$Length <- (as.numeric(Kidney4$Length))
Kidney4 <- Kidney4[!is.na(Kidney4$Length), ]
Kidney4 <- Kidney4[!is.na(Kidney4$value), ]
Kidney4 <- Kidney4[!is.infinite(Kidney4$Length), ]
Kidney4 <- Kidney4[!is.infinite(Kidney4$value), ]
Kidney4$Type <- "TTN"

# Subset Dystrophin
Kidney5 <- subset(Kidney, id=="ENSG00000198947")
Kidney5 <- melt(Kidney5, id.vars="Length", "foldChange")
Kidney5$value <- log2(as.numeric(Kidney5$value))
Kidney5$Length <- (as.numeric(Kidney5$Length))
Kidney5 <- Kidney5[!is.na(Kidney5$Length), ]
Kidney5 <- Kidney5[!is.na(Kidney5$value), ]
Kidney5 <- Kidney5[!is.infinite(Kidney5$Length), ]
Kidney5 <- Kidney5[!is.infinite(Kidney5$value), ]
Kidney5$Type <- "DMD"

# Merge and plot

Kidney6 <- rbind(Kidney2, Kidney3, Kidney4, Kidney5) 
Kidney_scatter <- ggplot(Kidney6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Kidney2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Kidney") + theme(legend.position="none")


# Process data for Liver

Liver <- Liver[grep("YES", Liver$coding), ]
Liver <- subset(Liver, baseMeanA>0, baseMeanB>0)
Liver <- subset(Liver, foldChange>0)
Liver2 <- melt(Liver, id.vars="Length", "foldChange")
Liver2$value <- log2(as.numeric(Liver2$value))
Liver2$Length <- (as.numeric(Liver2$Length))
Liver2 <- Liver2[!is.na(Liver2$Length), ]
Liver2 <- Liver2[!is.na(Liver2$value), ]
Liver2 <- Liver2[!is.infinite(Liver2$Length), ]
Liver2 <- Liver2[!is.infinite(Liver2$value), ]
Liver2$Type <- "Other"

# Subset recursive genes

Liver3 <- subset(Liver, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Liver3 <- melt(Liver3, id.vars="Length", "foldChange")
Liver3$value <- log2(as.numeric(Liver3$value))
Liver3$Length <- (as.numeric(Liver3$Length))
Liver3 <- Liver3[!is.na(Liver3$Length), ]
Liver3 <- Liver3[!is.na(Liver3$value), ]
Liver3 <- Liver3[!is.infinite(Liver3$Length), ]
Liver3 <- Liver3[!is.infinite(Liver3$value), ]
Liver3$Type <- "RSS"

# Subset Titin
Liver4 <- subset(Liver, id=="ENSG00000155657")
Liver4 <- melt(Liver4, id.vars="Length", "foldChange")
Liver4$value <- log2(as.numeric(Liver4$value))
Liver4$Length <- (as.numeric(Liver4$Length))
Liver4 <- Liver4[!is.na(Liver4$Length), ]
Liver4 <- Liver4[!is.na(Liver4$value), ]
Liver4 <- Liver4[!is.infinite(Liver4$Length), ]
Liver4 <- Liver4[!is.infinite(Liver4$value), ]
Liver4$Type <- "TTN"

# Subset Dystrophin
Liver5 <- subset(Liver, id=="ENSG00000198947")
Liver5 <- melt(Liver5, id.vars="Length", "foldChange")
Liver5$value <- log2(as.numeric(Liver5$value))
Liver5$Length <- (as.numeric(Liver5$Length))
Liver5 <- Liver5[!is.na(Liver5$Length), ]
Liver5 <- Liver5[!is.na(Liver5$value), ]
Liver5 <- Liver5[!is.infinite(Liver5$Length), ]
Liver5 <- Liver5[!is.infinite(Liver5$value), ]
Liver5$Type <- "DMD"

# Merge and plot

Liver6 <- rbind(Liver2, Liver3, Liver4, Liver5) 
Liver_scatter <- ggplot(Liver6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Liver2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Liver") + theme(legend.position="none")


# Process data for Lung

Lung <- Lung[grep("YES", Lung$coding), ]
Lung <- subset(Lung, baseMeanA>0, baseMeanB>0)
Lung <- subset(Lung, foldChange>0)
Lung2 <- melt(Lung, id.vars="Length", "foldChange")
Lung2$value <- log2(as.numeric(Lung2$value))
Lung2$Length <- (as.numeric(Lung2$Length))
Lung2 <- Lung2[!is.na(Lung2$Length), ]
Lung2 <- Lung2[!is.na(Lung2$value), ]
Lung2 <- Lung2[!is.infinite(Lung2$Length), ]
Lung2 <- Lung2[!is.infinite(Lung2$value), ]
Lung2$Type <- "Other"

# Subset recursive genes

Lung3 <- subset(Lung, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Lung3 <- melt(Lung3, id.vars="Length", "foldChange")
Lung3$value <- log2(as.numeric(Lung3$value))
Lung3$Length <- (as.numeric(Lung3$Length))
Lung3 <- Lung3[!is.na(Lung3$Length), ]
Lung3 <- Lung3[!is.na(Lung3$value), ]
Lung3 <- Lung3[!is.infinite(Lung3$Length), ]
Lung3 <- Lung3[!is.infinite(Lung3$value), ]
Lung3$Type <- "RSS"

# Subset Titin
Lung4 <- subset(Lung, id=="ENSG00000155657")
Lung4 <- melt(Lung4, id.vars="Length", "foldChange")
Lung4$value <- log2(as.numeric(Lung4$value))
Lung4$Length <- (as.numeric(Lung4$Length))
Lung4 <- Lung4[!is.na(Lung4$Length), ]
Lung4 <- Lung4[!is.na(Lung4$value), ]
Lung4 <- Lung4[!is.infinite(Lung4$Length), ]
Lung4 <- Lung4[!is.infinite(Lung4$value), ]
Lung4$Type <- "TTN"

# Subset Dystrophin
Lung5 <- subset(Lung, id=="ENSG00000198947")
Lung5 <- melt(Lung5, id.vars="Length", "foldChange")
Lung5$value <- log2(as.numeric(Lung5$value))
Lung5$Length <- (as.numeric(Lung5$Length))
Lung5 <- Lung5[!is.na(Lung5$Length), ]
Lung5 <- Lung5[!is.na(Lung5$value), ]
Lung5 <- Lung5[!is.infinite(Lung5$Length), ]
Lung5 <- Lung5[!is.infinite(Lung5$value), ]
Lung5$Type <- "DMD"

# Merge and plot

Lung6 <- rbind(Lung2, Lung3, Lung4, Lung5) 
Lung_scatter <- ggplot(Lung6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Lung2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Lung") + theme(legend.position="none")


# Process data for Lymph

Lymph <- Lymph[grep("YES", Lymph$coding), ]
Lymph <- subset(Lymph, baseMeanA>0, baseMeanB>0)
Lymph <- subset(Lymph, foldChange>0)
Lymph2 <- melt(Lymph, id.vars="Length", "foldChange")
Lymph2$value <- log2(as.numeric(Lymph2$value))
Lymph2$Length <- (as.numeric(Lymph2$Length))
Lymph2 <- Lymph2[!is.na(Lymph2$Length), ]
Lymph2 <- Lymph2[!is.na(Lymph2$value), ]
Lymph2 <- Lymph2[!is.infinite(Lymph2$Length), ]
Lymph2 <- Lymph2[!is.infinite(Lymph2$value), ]
Lymph2$Type <- "Other"

# Subset recursive genes

Lymph3 <- subset(Lymph, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Lymph3 <- melt(Lymph3, id.vars="Length", "foldChange")
Lymph3$value <- log2(as.numeric(Lymph3$value))
Lymph3$Length <- (as.numeric(Lymph3$Length))
Lymph3 <- Lymph3[!is.na(Lymph3$Length), ]
Lymph3 <- Lymph3[!is.na(Lymph3$value), ]
Lymph3 <- Lymph3[!is.infinite(Lymph3$Length), ]
Lymph3 <- Lymph3[!is.infinite(Lymph3$value), ]
Lymph3$Type <- "RSS"

# Subset Titin
Lymph4 <- subset(Lymph, id=="ENSG00000155657")
Lymph4 <- melt(Lymph4, id.vars="Length", "foldChange")
Lymph4$value <- log2(as.numeric(Lymph4$value))
Lymph4$Length <- (as.numeric(Lymph4$Length))
Lymph4 <- Lymph4[!is.na(Lymph4$Length), ]
Lymph4 <- Lymph4[!is.na(Lymph4$value), ]
Lymph4 <- Lymph4[!is.infinite(Lymph4$Length), ]
Lymph4 <- Lymph4[!is.infinite(Lymph4$value), ]
Lymph4$Type <- "TTN"

# Subset Dystrophin
Lymph5 <- subset(Lymph, id=="ENSG00000198947")
Lymph5 <- melt(Lymph5, id.vars="Length", "foldChange")
Lymph5$value <- log2(as.numeric(Lymph5$value))
Lymph5$Length <- (as.numeric(Lymph5$Length))
Lymph5 <- Lymph5[!is.na(Lymph5$Length), ]
Lymph5 <- Lymph5[!is.na(Lymph5$value), ]
Lymph5 <- Lymph5[!is.infinite(Lymph5$Length), ]
Lymph5 <- Lymph5[!is.infinite(Lymph5$value), ]
Lymph5$Type <- "DMD"

# Merge and plot

Lymph6 <- rbind(Lymph2, Lymph3, Lymph4, Lymph5) 
Lymph_scatter <- ggplot(Lymph6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Lymph2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Lymph") + theme(legend.position="none")


# Process data for Ovary

Ovary <- Ovary[grep("YES", Ovary$coding), ]
Ovary <- subset(Ovary, baseMeanA>0, baseMeanB>0)
Ovary <- subset(Ovary, foldChange>0)
Ovary2 <- melt(Ovary, id.vars="Length", "foldChange")
Ovary2$value <- log2(as.numeric(Ovary2$value))
Ovary2$Length <- (as.numeric(Ovary2$Length))
Ovary2 <- Ovary2[!is.na(Ovary2$Length), ]
Ovary2 <- Ovary2[!is.na(Ovary2$value), ]
Ovary2 <- Ovary2[!is.infinite(Ovary2$Length), ]
Ovary2 <- Ovary2[!is.infinite(Ovary2$value), ]
Ovary2$Type <- "Other"

# Subset recursive genes

Ovary3 <- subset(Ovary, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Ovary3 <- melt(Ovary3, id.vars="Length", "foldChange")
Ovary3$value <- log2(as.numeric(Ovary3$value))
Ovary3$Length <- (as.numeric(Ovary3$Length))
Ovary3 <- Ovary3[!is.na(Ovary3$Length), ]
Ovary3 <- Ovary3[!is.na(Ovary3$value), ]
Ovary3 <- Ovary3[!is.infinite(Ovary3$Length), ]
Ovary3 <- Ovary3[!is.infinite(Ovary3$value), ]
Ovary3$Type <- "RSS"

# Subset Titin
Ovary4 <- subset(Ovary, id=="ENSG00000155657")
Ovary4 <- melt(Ovary4, id.vars="Length", "foldChange")
Ovary4$value <- log2(as.numeric(Ovary4$value))
Ovary4$Length <- (as.numeric(Ovary4$Length))
Ovary4 <- Ovary4[!is.na(Ovary4$Length), ]
Ovary4 <- Ovary4[!is.na(Ovary4$value), ]
Ovary4 <- Ovary4[!is.infinite(Ovary4$Length), ]
Ovary4 <- Ovary4[!is.infinite(Ovary4$value), ]
Ovary4$Type <- "TTN"

# Subset Dystrophin
Ovary5 <- subset(Ovary, id=="ENSG00000198947")
Ovary5 <- melt(Ovary5, id.vars="Length", "foldChange")
Ovary5$value <- log2(as.numeric(Ovary5$value))
Ovary5$Length <- (as.numeric(Ovary5$Length))
Ovary5 <- Ovary5[!is.na(Ovary5$Length), ]
Ovary5 <- Ovary5[!is.na(Ovary5$value), ]
Ovary5 <- Ovary5[!is.infinite(Ovary5$Length), ]
Ovary5 <- Ovary5[!is.infinite(Ovary5$value), ]
Ovary5$Type <- "DMD"

# Merge and plot

Ovary6 <- rbind(Ovary2, Ovary3, Ovary4, Ovary5) 
Ovary_scatter <- ggplot(Ovary6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Ovary2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Ovary") + theme(legend.position="none")


# Process data for Prostate

Prostate <- Prostate[grep("YES", Prostate$coding), ]
Prostate <- subset(Prostate, baseMeanA>0, baseMeanB>0)
Prostate <- subset(Prostate, foldChange>0)
Prostate2 <- melt(Prostate, id.vars="Length", "foldChange")
Prostate2$value <- log2(as.numeric(Prostate2$value))
Prostate2$Length <- (as.numeric(Prostate2$Length))
Prostate2 <- Prostate2[!is.na(Prostate2$Length), ]
Prostate2 <- Prostate2[!is.na(Prostate2$value), ]
Prostate2 <- Prostate2[!is.infinite(Prostate2$Length), ]
Prostate2 <- Prostate2[!is.infinite(Prostate2$value), ]
Prostate2$Type <- "Other"

# Subset recursive genes

Prostate3 <- subset(Prostate, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Prostate3 <- melt(Prostate3, id.vars="Length", "foldChange")
Prostate3$value <- log2(as.numeric(Prostate3$value))
Prostate3$Length <- (as.numeric(Prostate3$Length))
Prostate3 <- Prostate3[!is.na(Prostate3$Length), ]
Prostate3 <- Prostate3[!is.na(Prostate3$value), ]
Prostate3 <- Prostate3[!is.infinite(Prostate3$Length), ]
Prostate3 <- Prostate3[!is.infinite(Prostate3$value), ]
Prostate3$Type <- "RSS"

# Subset Titin
Prostate4 <- subset(Prostate, id=="ENSG00000155657")
Prostate4 <- melt(Prostate4, id.vars="Length", "foldChange")
Prostate4$value <- log2(as.numeric(Prostate4$value))
Prostate4$Length <- (as.numeric(Prostate4$Length))
Prostate4 <- Prostate4[!is.na(Prostate4$Length), ]
Prostate4 <- Prostate4[!is.na(Prostate4$value), ]
Prostate4 <- Prostate4[!is.infinite(Prostate4$Length), ]
Prostate4 <- Prostate4[!is.infinite(Prostate4$value), ]
Prostate4$Type <- "TTN"

# Subset Dystrophin
Prostate5 <- subset(Prostate, id=="ENSG00000198947")
Prostate5 <- melt(Prostate5, id.vars="Length", "foldChange")
Prostate5$value <- log2(as.numeric(Prostate5$value))
Prostate5$Length <- (as.numeric(Prostate5$Length))
Prostate5 <- Prostate5[!is.na(Prostate5$Length), ]
Prostate5 <- Prostate5[!is.na(Prostate5$value), ]
Prostate5 <- Prostate5[!is.infinite(Prostate5$Length), ]
Prostate5 <- Prostate5[!is.infinite(Prostate5$value), ]
Prostate5$Type <- "DMD"

# Merge and plot

Prostate6 <- rbind(Prostate2, Prostate3, Prostate4, Prostate5) 
Prostate_scatter <- ggplot(Prostate6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Prostate2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Prostate") + theme(legend.position="none")


# Process data for Skeletal muscle

Skeletal <- Skeletal[grep("YES", Skeletal$coding), ]
Skeletal <- subset(Skeletal, baseMeanA>0, baseMeanB>0)
Skeletal <- subset(Skeletal, foldChange>0)
Skeletal2 <- melt(Skeletal, id.vars="Length", "foldChange")
Skeletal2$value <- log2(as.numeric(Skeletal2$value))
Skeletal2$Length <- (as.numeric(Skeletal2$Length))
Skeletal2 <- Skeletal2[!is.na(Skeletal2$Length), ]
Skeletal2 <- Skeletal2[!is.na(Skeletal2$value), ]
Skeletal2 <- Skeletal2[!is.infinite(Skeletal2$Length), ]
Skeletal2 <- Skeletal2[!is.infinite(Skeletal2$value), ]
Skeletal2$Type <- "Other"

# Subset recursive genes

Skeletal3 <- subset(Skeletal, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Skeletal3 <- melt(Skeletal3, id.vars="Length", "foldChange")
Skeletal3$value <- log2(as.numeric(Skeletal3$value))
Skeletal3$Length <- (as.numeric(Skeletal3$Length))
Skeletal3 <- Skeletal3[!is.na(Skeletal3$Length), ]
Skeletal3 <- Skeletal3[!is.na(Skeletal3$value), ]
Skeletal3 <- Skeletal3[!is.infinite(Skeletal3$Length), ]
Skeletal3 <- Skeletal3[!is.infinite(Skeletal3$value), ]
Skeletal3$Type <- "RSS"

# Subset Titin
Skeletal4 <- subset(Skeletal, id=="ENSG00000155657")
Skeletal4 <- melt(Skeletal4, id.vars="Length", "foldChange")
Skeletal4$value <- log2(as.numeric(Skeletal4$value))
Skeletal4$Length <- (as.numeric(Skeletal4$Length))
Skeletal4 <- Skeletal4[!is.na(Skeletal4$Length), ]
Skeletal4 <- Skeletal4[!is.na(Skeletal4$value), ]
Skeletal4 <- Skeletal4[!is.infinite(Skeletal4$Length), ]
Skeletal4 <- Skeletal4[!is.infinite(Skeletal4$value), ]
Skeletal4$Type <- "TTN"

# Subset Dystrophin
Skeletal5 <- subset(Skeletal, id=="ENSG00000198947")
Skeletal5 <- melt(Skeletal5, id.vars="Length", "foldChange")
Skeletal5$value <- log2(as.numeric(Skeletal5$value))
Skeletal5$Length <- (as.numeric(Skeletal5$Length))
Skeletal5 <- Skeletal5[!is.na(Skeletal5$Length), ]
Skeletal5 <- Skeletal5[!is.na(Skeletal5$value), ]
Skeletal5 <- Skeletal5[!is.infinite(Skeletal5$Length), ]
Skeletal5 <- Skeletal5[!is.infinite(Skeletal5$value), ]
Skeletal5$Type <- "DMD"

# Merge and plot

Skeletal6 <- rbind(Skeletal2, Skeletal3, Skeletal4, Skeletal5) 
Skeletal_scatter <- ggplot(Skeletal6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Skeletal2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Skeletal") + theme(legend.position="none")


# Process data for Testes

Testes <- Testes[grep("YES", Testes$coding), ]
Testes <- subset(Testes, baseMeanA>0, baseMeanB>0)
Testes <- subset(Testes, foldChange>0)
Testes2 <- melt(Testes, id.vars="Length", "foldChange")
Testes2$value <- log2(as.numeric(Testes2$value))
Testes2$Length <- (as.numeric(Testes2$Length))
Testes2 <- Testes2[!is.na(Testes2$Length), ]
Testes2 <- Testes2[!is.na(Testes2$value), ]
Testes2 <- Testes2[!is.infinite(Testes2$Length), ]
Testes2 <- Testes2[!is.infinite(Testes2$value), ]
Testes2$Type <- "Other"

# Subset recursive genes

Testes3 <- subset(Testes, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Testes3 <- melt(Testes3, id.vars="Length", "foldChange")
Testes3$value <- log2(as.numeric(Testes3$value))
Testes3$Length <- (as.numeric(Testes3$Length))
Testes3 <- Testes3[!is.na(Testes3$Length), ]
Testes3 <- Testes3[!is.na(Testes3$value), ]
Testes3 <- Testes3[!is.infinite(Testes3$Length), ]
Testes3 <- Testes3[!is.infinite(Testes3$value), ]
Testes3$Type <- "RSS"

# Subset Titin
Testes4 <- subset(Testes, id=="ENSG00000155657")
Testes4 <- melt(Testes4, id.vars="Length", "foldChange")
Testes4$value <- log2(as.numeric(Testes4$value))
Testes4$Length <- (as.numeric(Testes4$Length))
Testes4 <- Testes4[!is.na(Testes4$Length), ]
Testes4 <- Testes4[!is.na(Testes4$value), ]
Testes4 <- Testes4[!is.infinite(Testes4$Length), ]
Testes4 <- Testes4[!is.infinite(Testes4$value), ]
Testes4$Type <- "TTN"

# Subset Dystrophin
Testes5 <- subset(Testes, id=="ENSG00000198947")
Testes5 <- melt(Testes5, id.vars="Length", "foldChange")
Testes5$value <- log2(as.numeric(Testes5$value))
Testes5$Length <- (as.numeric(Testes5$Length))
Testes5 <- Testes5[!is.na(Testes5$Length), ]
Testes5 <- Testes5[!is.na(Testes5$value), ]
Testes5 <- Testes5[!is.infinite(Testes5$Length), ]
Testes5 <- Testes5[!is.infinite(Testes5$value), ]
Testes5$Type <- "DMD"

# Merge and plot

Testes6 <- rbind(Testes2, Testes3, Testes4, Testes5) 
Testes_scatter <- ggplot(Testes6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Testes2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Testes") + theme(legend.position="none")


# Process data for Thyroid

Thyroid <- Thyroid[grep("YES", Thyroid$coding), ]
Thyroid <- subset(Thyroid, baseMeanA>0, baseMeanB>0)
Thyroid <- subset(Thyroid, foldChange>0)
Thyroid2 <- melt(Thyroid, id.vars="Length", "foldChange")
Thyroid2$value <- log2(as.numeric(Thyroid2$value))
Thyroid2$Length <- (as.numeric(Thyroid2$Length))
Thyroid2 <- Thyroid2[!is.na(Thyroid2$Length), ]
Thyroid2 <- Thyroid2[!is.na(Thyroid2$value), ]
Thyroid2 <- Thyroid2[!is.infinite(Thyroid2$Length), ]
Thyroid2 <- Thyroid2[!is.infinite(Thyroid2$value), ]
Thyroid2$Type <- "Other"

# Subset recursive genes

Thyroid3 <- subset(Thyroid, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
Thyroid3 <- melt(Thyroid3, id.vars="Length", "foldChange")
Thyroid3$value <- log2(as.numeric(Thyroid3$value))
Thyroid3$Length <- (as.numeric(Thyroid3$Length))
Thyroid3 <- Thyroid3[!is.na(Thyroid3$Length), ]
Thyroid3 <- Thyroid3[!is.na(Thyroid3$value), ]
Thyroid3 <- Thyroid3[!is.infinite(Thyroid3$Length), ]
Thyroid3 <- Thyroid3[!is.infinite(Thyroid3$value), ]
Thyroid3$Type <- "RSS"

# Subset Titin
Thyroid4 <- subset(Thyroid, id=="ENSG00000155657")
Thyroid4 <- melt(Thyroid4, id.vars="Length", "foldChange")
Thyroid4$value <- log2(as.numeric(Thyroid4$value))
Thyroid4$Length <- (as.numeric(Thyroid4$Length))
Thyroid4 <- Thyroid4[!is.na(Thyroid4$Length), ]
Thyroid4 <- Thyroid4[!is.na(Thyroid4$value), ]
Thyroid4 <- Thyroid4[!is.infinite(Thyroid4$Length), ]
Thyroid4 <- Thyroid4[!is.infinite(Thyroid4$value), ]
Thyroid4$Type <- "TTN"

# Subset Dystrophin
Thyroid5 <- subset(Thyroid, id=="ENSG00000198947")
Thyroid5 <- melt(Thyroid5, id.vars="Length", "foldChange")
Thyroid5$value <- log2(as.numeric(Thyroid5$value))
Thyroid5$Length <- (as.numeric(Thyroid5$Length))
Thyroid5 <- Thyroid5[!is.na(Thyroid5$Length), ]
Thyroid5 <- Thyroid5[!is.na(Thyroid5$value), ]
Thyroid5 <- Thyroid5[!is.infinite(Thyroid5$Length), ]
Thyroid5 <- Thyroid5[!is.infinite(Thyroid5$value), ]
Thyroid5$Type <- "DMD"

# Merge and plot

Thyroid6 <- rbind(Thyroid2, Thyroid3, Thyroid4, Thyroid5) 
Thyroid_scatter <- ggplot(Thyroid6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(Thyroid2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. Thyroid") + theme(legend.position="none")


# Process data for White blood cells

White <- White[grep("YES", White$coding), ]
White <- subset(White, baseMeanA>0, baseMeanB>0)
White <- subset(White, foldChange>0)
White2 <- melt(White, id.vars="Length", "foldChange")
White2$value <- log2(as.numeric(White2$value))
White2$Length <- (as.numeric(White2$Length))
White2 <- White2[!is.na(White2$Length), ]
White2 <- White2[!is.na(White2$value), ]
White2 <- White2[!is.infinite(White2$Length), ]
White2 <- White2[!is.infinite(White2$value), ]
White2$Type <- "Other"

# Subset recursive genes

White3 <- subset(White, id=="ENSG00000182985" | id=="ENSG00000175161" | id=="ENSG00000151150" | id=="ENSG00000113448" | id=="ENSG00000183715" | id=="ENSG00000185008" | id=="ENSG00000185352" | id =="ENSG00000149294")
White3 <- melt(White3, id.vars="Length", "foldChange")
White3$value <- log2(as.numeric(White3$value))
White3$Length <- (as.numeric(White3$Length))
White3 <- White3[!is.na(White3$Length), ]
White3 <- White3[!is.na(White3$value), ]
White3 <- White3[!is.infinite(White3$Length), ]
White3 <- White3[!is.infinite(White3$value), ]
White3$Type <- "RSS"

# Subset Titin
White4 <- subset(White, id=="ENSG00000155657")
White4 <- melt(White4, id.vars="Length", "foldChange")
White4$value <- log2(as.numeric(White4$value))
White4$Length <- (as.numeric(White4$Length))
White4 <- White4[!is.na(White4$Length), ]
White4 <- White4[!is.na(White4$value), ]
White4 <- White4[!is.infinite(White4$Length), ]
White4 <- White4[!is.infinite(White4$value), ]
White4$Type <- "TTN"

# Subset Dystrophin
White5 <- subset(White, id=="ENSG00000198947")
White5 <- melt(White5, id.vars="Length", "foldChange")
White5$value <- log2(as.numeric(White5$value))
White5$Length <- (as.numeric(White5$Length))
White5 <- White5[!is.na(White5$Length), ]
White5 <- White5[!is.na(White5$value), ]
White5 <- White5[!is.infinite(White5$Length), ]
White5 <- White5[!is.infinite(White5$value), ]
White5$Type <- "DMD"

# Merge and plot

White6 <- rbind(White2, White3, White4, White5) 
White_scatter <- ggplot(White6, aes(Length, value, color=Type, size=factor(Type), fill=factor(Type))) + geom_point() + scale_colour_manual(values=c("black", "grey65", "red", "blue")) + scale_size_manual(values=c(1.6, 0.8, 1.6, 1.6)) + scale_x_log10() + geom_hline(yintercept=0, colour="steelblue", linetype = "longdash", size=0.5) + geom_vline(xintercept=150000, colour="steelblue", linetype = "longdash", size=0.5) + stat_smooth(data=subset(White2, Length>0), size=0.7, colour="green3") + ylim(c(-20,20)) + labs(x="log10(Gene Length)", y="log2(foldChange)", title="Brain vs. White") + theme(legend.position="none")



pdf(file='fig/ED1A_Bodymap_vs_brain_GeneLength.pdf', width=9, height=11)
grid.arrange(Adipose_scatter, Adrenal_scatter, Breast_scatter, Colon_scatter, Heart_scatter, Kidney_scatter, Liver_scatter, Lung_scatter, Lymph_scatter, Ovary_scatter, Prostate_scatter, Skeletal_scatter, Testes_scatter, Thyroid_scatter, White_scatter, ncol=3)
dev.off()


