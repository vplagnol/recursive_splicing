source('scripts/tools.R')
cassette <- TRUE

data <- read.csv('data/ensembl_annotations_complete.csv', stringsAsFactors = FALSE)
data$ID <- paste('J', data$ID, sep = '')

my.bams <- scan('data/BAMlist.tab', what = 'character')
pos <- NULL


                       
RPKM.table <- read.table('data/RPKM_brain_processed.tab', header = TRUE, sep = '\t')


bed.file <- 'data/group_2771_bestof-FUS-C35_sum_G_hg19--ensembl59_from_439-685-686_bedGraph-cDNA-hits-in-genome.bed'
bed.file <- read.table(file = bed.file, sep = '\t', skip = 1, header = FALSE)


#### get a list of introns
data$intronID <- paste(data$Chromosome, '_', data$IntronStart, '_', data$IntronEnd, sep = '')
list.introns <- subset( data[, c('Chromosome', 'IntronStart', 'IntronEnd', 'intronID', 'Gene', 'HUGO', 'Strand')], data$dist.to.exons <= 4)
list.introns <- do.call(rbind.data.frame, by (list.introns, INDICES = list.introns$intronID, FUN = function(x) {x[1,]}))
list.introns$RPKM <- RPKM.table$RPKM [ match(list.introns$Gene, table = RPKM.table$ensemblID) ]

############ can we look at cassette exons?
if (cassette) {
  list.introns <- read.csv('data/cassette_exons_with_long_introns.csv', stringsAsFactors = FALSE)
  list.introns$Gene <- list.introns$intronRightGene
  list.introns$HUGO <- list.introns$intronRightHUGO
  list.introns$IntronStart <- pmin( list.introns$intronStartRight, list.introns$intronEndRight)
  list.introns$IntronEnd <- pmax( list.introns$intronStartRight, list.introns$intronEndRight)
  list.introns$intronID <- list.introns$intronIDRight
  list.introns$Chromosome <- gsub(pattern = '^chr', replacement = '', list.introns$X.chrom)
  list.introns$Strand <- ifelse ( list.introns$strand == '-', -1, +1)
  list.introns$ID <- paste('ID', 1:nrow(list.introns), sep = '')
}




##############

if (!cassette) {
  choice.list <- 1
  if (choice.list == 1) junction.set <- 'data/Regression_IDs_1_Exon_Intron.csv'
  if (choice.list == 2) junction.set <- 'data/Regression_IDs_2_Intron_Exon.csv'
  list.good.junctions <- scan(junction.set, what = 'character')
  output.table.FUS <- paste('data/selected_junctions_with_exons_FUS_list', choice.list, '.csv', sep = '')
  output.table.RNAseq <- paste('data/selected_junctions_with_exons_RNAseq_list', choice.list, '.csv', sep = '')
}

if (cassette) {
  choice.list <- 'cassette'
  output.table.FUS <- paste('data/selected_junctions_with_exons_FUS_list', choice.list, '.csv', sep = '')
  output.table.RNAseq <- paste('data/selected_junctions_with_exons_RNAseq_list', choice.list, '.csv', sep = '')
}




#list.introns <- subset(list.introns, HUGO == 'NTM')

#list.introns <- subset(list.introns, Gene %in% RPKM.table$ensemblID)  ##we need a RPKM value to get the slopes right
#list.introns <- subset(list.introns, HUGO %in% c('HPSE2'))[1,]
#list.introns <- subset(list.introns, HUGO == 'MAGI1')
#list.introns <- subset(list.introns, HUGO == 'OPCML')

nintrons <- nrow(list.introns)
list.introns$RPKM <- NA
list.introns.with.exons <- list.introns
list.introns.with.exons.FUS <- list.introns
list.introns.without.exons <- list.introns


#known.exons <- read.csv('data/all_exon_intron_junctions.tab', header = FALSE, col.names = c('chromosome', 'position'), sep = '\t')  ##ensembl exons
known.exons <- read.table('data/new_1_exons_genes.sort.uniq.bed', col.names = c('chromosome', 'start', 'end', 'gene'), sep = '\t', stringsAsFactors = FALSE)  ##UCSC exons I think
known.exons$chromosome <- gsub(pattern = 'chr', replacement = '', known.exons$chromosome)
known.exons.1 <- data.frame(chromosome  = known.exons$chromosome, position = known.exons$start)
known.exons.2 <- data.frame(chromosome  = known.exons$chromosome, position = known.exons$end)
known.exons <- rbind.data.frame( known.exons.1, known.exons.2 )



################ Now look at each intron
final.RNAseq <- data.frame()
final.FUS <- data.frame()

for (i in 1:nintrons) {
  
  RPKM.value <- subset(RPKM.table, ensemblID == list.introns$Gene[i], 'RPKM', drop = TRUE)
  if (length(RPKM.value) == 0) RPKM.value <- 1
  
  
  list.introns$RPKM[ i ] <- RPKM.value
  message(list.introns$HUGO[i], ' ', i, ' ', nintrons, ' ', RPKM.value)



  if (cassette) output.pdf <- paste('fig/seesaw/cassette_', list.introns$HUGO[ i ], '_', list.introns$intronID[ i], '_list', choice.list, '.pdf', sep = '') else {
    output.pdf <- paste('fig/seesaw/', list.introns$HUGO[ i ], '_', list.introns$intronID[ i], '_list', choice.list, '.pdf', sep = '')
  }


##########
  output.reads <- paste('temp/reads/', list.introns$intronID[ i], sep = '')
  region <- paste('chr', list.introns$Chromosome[i], ':', list.introns$IntronStart[ i] , '-', list.introns$IntronEnd[ i ], sep = '')
  if (!file.exists(output.reads)) {
    for (bam in my.bams) {
      my.cmd <- paste('/cluster/project8/vyp/vincent/Software/samtools-0.1.18/samtools view', bam, region, " | awk '{if (( $9 > 0) && ($9 < 400)) print ($4 + $8)}' >> ", output.reads)
      system(my.cmd)
    }
  }
  pos.RNAseq <- as.numeric(scan(output.reads, what = 'numeric'))
  pos.RNAseq <- as.numeric(pos.RNAseq)/2000
  
##########
  bed.file.loc <- subset(bed.file, V1 == paste('chr', list.introns$Chromosome[ i ], sep = ''))
  bed.file.loc <- subset(bed.file.loc, V2 > list.introns$IntronStart[ i ] & V3 < list.introns$IntronEnd[ i ])
  if (list.introns$Strand[i] == 1) bed.file.loc <- subset(bed.file.loc, V4 > 0)
  if (list.introns$Strand[i] == -1) bed.file.loc <- subset(bed.file.loc, V4 < 0)
  
  pos.FUS <-  rep(bed.file.loc$V2 + bed.file.loc$V3, times = abs(bed.file.loc$V4))/2000


######
  length.kb <- floor((list.introns$IntronEnd[ i] - list.introns$IntronStart[i])/1000)
  if (cassette) {
    junctions <- data.frame( ID = list.introns$ID[ i ],  warren.position.junction = list.introns$junctionPos[ i ] , stringsAsFactors = FALSE)
  } else {
    junctions <- subset(data, warren.position.junction >= list.introns$IntronStart[ i ] &
                        warren.position.junction <= list.introns$IntronEnd[ i ]
                        & Chromosome == list.introns$Chromosome[ i ] & dist.to.exons <= 4)
    
    if (choice.list == 1) junctions <- subset(junctions, abs(warren.position.exon - list.introns$IntronStart[ i ] < 10) | abs(warren.position.exon - list.introns$IntronEnd[ i ] < 10) )
    junctions <- subset(junctions, ID %in% list.good.junctions)
  }    
  message('Number of junctions ', nrow(junctions))
  goforit.FUS <- (length(pos.FUS) > 100) && (nrow(junctions) > 0)
  goforit.RNAseq <- (length(pos.RNAseq) > 100) && (nrow(junctions) > 0 )

  if (goforit.RNAseq) {
    unique.junctions <- do.call(rbind.data.frame, by( junctions[, c('ID', 'warren.position.junction')], INDICES = junctions$ID, function(x) {x[1,]}))

    #list.introns.with.exons$pdf[ i ] <- output.pdf
  
    my.xlim <- c(list.introns$IntronStart[ i ]/1000, (list.introns$IntronEnd[ i ])/1000)
    length.kb <- my.xlim[2] - my.xlim[1]
    pos.RNAseq <- subset(pos.RNAseq, pos.RNAseq > my.xlim[1] & pos.RNAseq < my.xlim[2])
    
    if (cassette) {
      known.exons.loc <- c()
    } else {
      known.exons.loc <- subset(known.exons, chromosome == list.introns$Chromosome[ i ] & position/1000 > my.xlim[1] & position/1000 < my.xlim[2], 'position', drop = TRUE)/1000
    }
    
    njunc.loc <- nrow(unique.junctions)
    over.P <- c()
    for (z in 1:njunc.loc) {
      res.with <- get.slope.estimate (pos.reads.kb = pos.RNAseq,
                                      start.intron = my.xlim[1],
                                      end.intron = my.xlim[2],
                                      strand = list.introns$Strand[ i ],
                                      warren.junctions.pos = unique.junctions$warren.position.junction[ z ] / 1000,
                                      warren.junctions.IDs = unique.junctions$ID[z],
                                      extra.exons.pos = known.exons.loc,
                                      my.title = paste(output.pdf, ', with additional exons included'),
                                      plot = FALSE)
      res.with$pdf <- output.pdf
      list.introns.with.exons <- copy.list ( list = res.with, data.frame = list.introns.with.exons, row = i)
      final.RNAseq <- rbind.data.frame(final.RNAseq, list.introns.with.exons[i,])
      over.P <- c(over.P, res.with$my.overall.P)

      message('And now FUS')
      res.with <- get.slope.estimate (pos.reads.kb = pos.FUS,
                                      start.intron = my.xlim[1],
                                      end.intron = my.xlim[2],
                                      strand = list.introns$Strand[ i ],
                                      warren.junctions.pos = unique.junctions$warren.position.junction[ z ] / 1000,
                                      warren.junctions.IDs = unique.junctions$ID[z],
                                      extra.exons.pos = known.exons.loc,
                                      my.title = paste('FUS, ', basename(output.pdf), ', with additional exons included'),
                                      plot = FALSE)
      res.with$pdf <- output.pdf
      list.introns.with.exons.FUS <- copy.list ( list = res.with, data.frame = list.introns.with.exons.FUS, row = i)
      final.FUS <- rbind.data.frame(final.FUS, list.introns.with.exons.FUS[i,])
    }


############
    z <- which.min(over.P)
    pdf(output.pdf, width = 9, height= 9)
    par(mar = c(6, 4, 4, 1), mfrow = c(2, 1))

    res.with <- get.slope.estimate (pos.reads.kb = pos.RNAseq,
                                    start.intron = my.xlim[1],
                                    end.intron = my.xlim[2],
                                    strand = list.introns$Strand[ i ],
                                    warren.junctions.pos = unique.junctions$warren.position.junction[ z ] / 1000,
                                    warren.junctions.IDs = unique.junctions$ID[z],
                                    extra.exons.pos = known.exons.loc,
                                    my.title = paste('RNA-Seq, ', basename(output.pdf), ', with additional exons included'),
                                    plot = TRUE)
    list.introns.with.exons <- copy.list ( list = res.with, data.frame = list.introns.with.exons, row = i)
    
    message('And now plort for FUS')
    res.with <- get.slope.estimate (pos.reads.kb = pos.FUS,
                                    start.intron = my.xlim[1],
                                    end.intron = my.xlim[2],
                                    strand = list.introns$Strand[ i ],
                                    warren.junctions.pos = unique.junctions$warren.position.junction[ z ] / 1000,
                                    warren.junctions.IDs = unique.junctions$ID[z],
                                    extra.exons.pos = known.exons.loc,
                                    my.title = paste('FUS, ', basename(output.pdf), ', with additional exons included'),
                                    plot = TRUE)
    dev.off()
    print(output.pdf)

    
    write.csv(x = final.RNAseq, row.names = FALSE, file = output.table.RNAseq) ##print the output after each gene
    write.csv(x = final.FUS, row.names = FALSE, file = output.table.FUS) ##print the output after each gene
  }

}




