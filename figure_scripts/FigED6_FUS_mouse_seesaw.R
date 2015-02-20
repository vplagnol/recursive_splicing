##creates histograms with regression lines for mouse FUS data- Extended Data 6

library(GenomicRanges)

count.data <- read.table("data/group_3094_FUS-E18-all_sum_G_mm9--ensembl59_from_2320-3089-3090_bedGraph-cDNA-hits-in-genome.bed.gz", header = FALSE, sep = "\t", skip = 1)

intron.positions <- read.csv("data/Mouse_RSS_introns.csv", header = FALSE, stringsAsFactors = FALSE)
mouse.RSS <- read.csv("data/Mouse_RSS_mm9.csv", header = FALSE, stringsAsFactors = FALSE)


output.pdf <- "fig/ED6_mouse_FUS_seesaw.pdf"
pdf(output.pdf, width = 4*2, height = 3*4)
par(mfrow = c(4, 2))
n.sites <- nrow(mouse.RSS)
region.halfwidth <- 200*1000

pval.tables <- data.frame(clean.names = rep(NA_character_, 8), ANOVA.pvals = NA, stringsAsFactors = FALSE)

for (i in 1:8) {
  message("Plotting RS ", i)
  chromosome <- mouse.RSS$V2[ i ]
  left <- mouse.RSS$V3[ i ] - region.halfwidth
  right <- mouse.RSS$V3[ i ] + region.halfwidth
  RSS.name <- gsub(pattern = "_1$", replacement = "", mouse.RSS$V1[ i ])

  left <- intron.positions$V3[ match(tolower(RSS.name), tolower(intron.positions$V1)) ]
  right <- intron.positions$V4[ match(tolower(RSS.name), tolower(intron.positions$V1)) ]
  clean.name <- intron.positions$V1[ match(tolower(RSS.name), tolower(intron.positions$V1)) ]
  
  count.data.loc <- subset(count.data, V1 == chromosome & V2 >= left & V2 <= right)
  G.count.data.loc <- GRanges (seqnames = count.data.loc$V1, IRanges( count.data.loc$V2, count.data.loc$V3)) 

  my.grid <- GRanges( seqnames = chromosome, IRanges( start =  seq(from = left, to = right - 5000, by = 5000), end =  seq(from = left, to = right - 5000, by = 5000) + 5000))

  overlaps <- findOverlaps (query = G.count.data.loc, subject = my.grid)

  combined.sum <- tapply(abs(count.data.loc$V4[ overlaps@queryHits ]),
                         INDEX =  factor( overlaps@subjectHits, levels = 1:length(my.grid)),
                         FUN = sum)
  
  for.regression <- data.frame(x = 0.5*(start(my.grid) + end(my.grid)), y = ifelse ( is.na(combined.sum), 0, combined.sum) )  ##check the NA issue, looks suspicious
  for.regression$RSS <- for.regression$x > mouse.RSS$V3[ i ] 

  if (RSS.name == "PDE4D") for.regression$sec.RSS <- for.regression$x > mouse.RSS$V3[ 9 ]
  if (RSS.name == "CADM2") for.regression$sec.RSS <- for.regression$x > 67143137  ### adds the second RSS site in CADM2
  
  artificial.reads <- rep(x = for.regression$x / 1000, times = for.regression$y) ## this is a bit of a hack but well... I want a histogram
  my.hist <- hist( artificial.reads, breaks = seq(from = left, to = right , by = 5000) / 1000, plot = FALSE)
  
  for.regression$xclean <- my.hist$mids - min(my.hist$mids)
  
  ####### exclude the region outside of genes- no need anymore because I use proper intron definitions
  #for.regression$exclude <- cumsum(for.regression$y) - max(cumsum(for.regression$y)) >= -10 | cumsum(for.regression$y) < 10
  #for.regression <- subset(for.regression, ! exclude )
  
  basic.model <- lm (for.regression, formula = " y ~ xclean")
  with.RSS <- lm (for.regression, formula = " y ~ xclean + RSS")

  print(anova(basic.model, with.RSS))
  pval.tables$clean.names[i]  <- clean.name
  pval.tables$ANOVA.pvals[i]  <- anova(basic.model, with.RSS)[2,6]
  
  if (RSS.name %in% c("PDE4D", "CADM2")) {  with.RSS <- lm (for.regression, formula = " y ~ xclean + RSS + sec.RSS" )}
  
  fitted.ordered <- fitted(basic.model)
  fitted.with.RSS <- fitted( with.RSS )

  label <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
  
  plot(my.hist,
       xlab = 'Position (kb)',
       ylab = 'Read count',
       xlim = c(left, right)/1000,
       col = 'red',
       main = paste0(label[ i ], ". ", clean.name))
  
  points ( x = for.regression$x/ 1000, y = fitted.ordered, type = 'l', lty = 1, col = 'black', lwd = 2)
  points ( x = for.regression$x/ 1000, y = fitted.with.RSS, type = 'l', lty = 1, col = 'blue', lwd = 2)
  abline ( v = mouse.RSS$V3[ i ]/ 1000, col = "blue")


  if (RSS.name == "PDE4D") {  abline ( v = mouse.RSS$V3[ 9 ]/ 1000, col = "blue")}
  if (RSS.name == "CADM2") {  abline ( v = 67143137/1000, col = "blue")}
}


dev.off()
print(output.pdf)
