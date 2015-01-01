library(GenomicRanges)
library(biomaRt)
source('seesaw_pattern_detection/tools.R')

input.warren <- 'data/Final_150K_junction_analysis.csv'
#input.warren <-'data/fromWarren/Whole_genome_50Kplus_updated_15_07_2013_v1_rec_details_v08_trunc.txt'

  
#system('cut -f1-14 data/fromWarren/Whole_genome_50Kplus_updated_15_07_2013_v1_rec_details.txt > temp/warren.tab')
Warren.table <- read.csv(input.warren, stringsAsFactors = FALSE)
ensembl.table <- read.table('data/transcript_info.tab', header = TRUE, na.string = '-9')

### remove NAs for now so that the GRanges can work
ensembl.table$IntronStart <- ifelse ( ensembl.table$nIntrons == 0, -1, ensembl.table$IntronStart)
ensembl.table$IntronEnd <- ifelse ( ensembl.table$nIntrons == 0, 0, ensembl.table$IntronEnd)

### Now look for overlap
introns.GRanges <- GRanges(seqnames = ensembl.table$Chromosome,
                           IRanges(start= ensembl.table$IntronStart, end = ensembl.table$IntronEnd))


first.intron <- Warren.table$First.part.type == 'Intron'
Warren.table$warren.position.junction <- ifelse ( first.intron, Warren.table$First.part.coordinate, Warren.table$Second.part.coordinate)  ##that is the junction position
Warren.table$warren.position.exon <- ifelse ( !first.intron, Warren.table$First.part.coordinate, Warren.table$Second.part.coordinate)    ##that is the exon position


warren.GRanges.introns <- GRanges(seqnames = gsub(pattern = 'chr', replacement = '', Warren.table$Chromosome),
                                  IRanges(start = Warren.table$warren.position.junction, end = Warren.table$warren.position.junction + 1))


#### check which are in UCSC or ensembl junctions
UCSC.exons <- read.table('data/new_1_exons_genes.sort.uniq.bed', col.names = c('chromosome', 'start', 'end', 'gene'), sep = '\t', stringsAsFactors = FALSE)
UCSC.exons$chromosome <- gsub(pattern = 'chr', replacement = '', UCSC.exons$chromosome)
UCSC.exons.1 <- data.frame(chromosome  = UCSC.exons$chromosome, position = UCSC.exons$start)
UCSC.exons.2 <- data.frame(chromosome  = UCSC.exons$chromosome, position = UCSC.exons$end)
UCSC.exons <- rbind.data.frame( UCSC.exons.1, UCSC.exons.2 )

UCSC.GRanges.introns <- GRanges(seqnames = UCSC.exons$chromosome,
                                  IRanges(start = UCSC.exons$position -4 , end = UCSC.exons$position + 4))
test1 <- findOverlaps(query = warren.GRanges.introns, subject = UCSC.GRanges.introns) ##look for overlap between intron and Warren's junctions
Warren.table$junction.near.UCSC.exon.junction <- FALSE
Warren.table$junction.near.UCSC.exon.junction  [test1@queryHits ] <- TRUE

### and now ensembl
ensembl.exons <- read.csv('data/all_ensembl_exon_intron_junctions.tab', header = FALSE, col.names = c('chromosome', 'position'), sep = '\t')  ##ensembl exons
ensembl.GRanges.introns <- GRanges(seqnames = ensembl.exons$chromosome,
                                IRanges(start = ensembl.exons$position -4 , end = ensembl.exons$position + 4))
test2 <- findOverlaps(query = warren.GRanges.introns, subject = ensembl.GRanges.introns) ##look for overlap between intron and Warren's junctions
Warren.table$junction.near.ensembl.exon.junction <- FALSE
Warren.table$junction.near.ensembl.exon.junction  [test2@queryHits ] <- TRUE


##look for overlap between intron and Warren's junctions
test <- findOverlaps(query = introns.GRanges, subject = warren.GRanges.introns)
test <- data.frame( queryHits = test@queryHits, subjectHits = test@subjectHits)
test.v2 <- data.frame ( queryHits = subset(1:nrow(ensembl.table), !1:nrow(ensembl.table) %in% test$queryHits), subjectHits = NA)
test <- rbind.data.frame (test, test.v2)
test <- test[ order(test$queryHits), ]

####

keep.warren <- subset(names(Warren.table), ! names(Warren.table) %in% c('Chromosome', 'First.part.coordinate', 'Second.part.coordinate' ,'Gene.name', 'Entrez.ID', 'First.part.type', 'Second.part.type', 'Strand'))
Warren.table.v2 <- data.frame (  Chromosome = Warren.table[ test$subjectHits , 'Chromosome'])
for (i in names(ensembl.table)) {Warren.table.v2[, i ] <- ensembl.table[ test$queryHits , i]}
for (i in keep.warren) {Warren.table.v2[, i ] <- Warren.table[ test$subjectHits , i]}
Warren.table.v2$RS.site <- !is.na(Warren.table.v2$Junction.reported.within.418nt.)
Warren.table.v2$dist.to.exons <- pmin( abs(Warren.table.v2$warren.position.exon - Warren.table.v2$IntronStart), abs(Warren.table.v2$warren.position.exon - Warren.table.v2$IntronEnd))
Warren.table.v2 <- subset( Warren.table.v2, dist.to.exons < 4) ### now check whether we are indeed next to the end of that intron

ensembl.table <- Warren.table.v2
ensembl.table$IntronLength.afterRS <- ifelse (ensembl.table$RS.site,
                                              pmax( abs(ensembl.table$IntronStart - ensembl.table$warren.position.junction), abs(ensembl.table$IntronEnd - ensembl.table$warren.position.junction)),
                                              ensembl.table$IntronLength)

#### put the NAs back where they should be
ensembl.table$IntronStart <- ifelse ( ensembl.table$nIntrons == 0, NA, ensembl.table$IntronStart)
ensembl.table$IntronEnd <- ifelse ( ensembl.table$nIntrons == 0, NA, ensembl.table$IntronEnd)



########## Find the matching mouse gene

species <- 'hsapiens_gene_ensembl'
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
my.db1 <- getBM(attributes= c('ensembl_gene_id', 'mmusculus_homolog_ensembl_gene'), mart=ensembl)
ensembl.table$Mouse.ensemblID <- my.db1$mmusculus_homolog_ensembl_gene[ match( table = my.db1$ensembl_gene_id, ensembl.table$Gene) ]
#attributes <- listAttributes(ensembl)

### Find a proper name for the mouse gene
species <- 'mmusculus_gene_ensembl'
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
my.db2 <- getBM(attributes= c('ensembl_gene_id', 'external_gene_id'), mart=ensembl)
ensembl.table$Mouse.HUGO <- my.db2$external_gene_id[ match( table = my.db2$ensembl_gene_id, ensembl.table$Mouse.ensemblID) ]


### cell type information
cell.type <- read.table('data/cell_type.txt', sep = '\t', header = TRUE)
ensembl.table$cell.type.from.Jernej <- cell.type$cell.type[ match( ensembl.table$Mouse.HUGO, table = cell.type$Symbol) ]

### Add the FUS/TDP43 data
FUS <- read.table('/cluster/project9/vyp2/Cleveland_TDP43_FUS_RNASeq/Cleveland/FUS/processed/deseq/control_KD/deseq_FUS_differential_expression.tab', header = TRUE)
FUS <- add.zscore(FUS, set = 'FUS')
names(FUS)[1] <- 'Mouse.ensemblID'
names(FUS) <- replace(names(FUS), names(FUS) == 'zscore', 'zscore.FUS.Cleveland')
names(FUS) <- replace(names(FUS), names(FUS) == 'conditionKD', 'conditionKD.FUS')
ensembl.table <- merge(ensembl.table, FUS[, c('Mouse.ensemblID', 'zscore.FUS.Cleveland', 'conditionKD.FUS') ], by = 'Mouse.ensemblID', all.x = TRUE)

TDP43 <- read.table('/cluster/project9/vyp2/Cleveland_TDP43_FUS_RNASeq/Cleveland/TDP43/processed/deseq/control_KD/deseq_TDP43_differential_expression.tab', header = TRUE)
TDP43 <- add.zscore(TDP43, set = 'TDP43')
names(TDP43)[1] <- 'Mouse.ensemblID'
names(TDP43) <- replace(names(TDP43), names(TDP43) == 'zscore', 'zscore.TDP43.Cleveland')
names(TDP43) <- replace(names(TDP43), names(TDP43) == 'conditionKD', 'conditionKD.TDP43')
ensembl.table <- merge(ensembl.table, TDP43[, c('Mouse.ensemblID', 'zscore.TDP43.Cleveland', 'conditionKD.TDP43') ], by = 'Mouse.ensemblID', all.x = TRUE)

write.csv(x = ensembl.table, file = 'data/ensembl_annotations_complete.csv', row.names = FALSE)


#### Some basic stats
print(table(subset(ensembl.table, IntronLength > 150000)$RS.site))
print(table(subset(ensembl.table, zscore.FUS.Cleveland < -3 & IntronLength > 150000)$RS.site))
print(table(subset(ensembl.table, zscore.TDP43.Cleveland < -3 & IntronLength > 150000)$RS.site))
