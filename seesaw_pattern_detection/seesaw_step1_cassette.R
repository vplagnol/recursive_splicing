library(GenomicRanges)
library(biomaRt)
source('scripts/tools.R')

minlength <- 25000

cassette.exons <- read.csv('data/CassetteExons.csv', stringsAsFactors = FALSE)
ensembl.table <- read.table('data/transcript_info.tab', header = TRUE, na.string = '-9')
ensembl.table <- subset(ensembl.table, !is.na(IntronEnd))
ensembl.table <- subset(ensembl.table, !is.na(IntronStart))
ensembl.table <- subset(ensembl.table, IntronLength > minlength)
ensembl.table$cleanStart <- pmin(ensembl.table$IntronStart,ensembl.table$IntronEnd)
ensembl.table$cleanEnd <- pmax(ensembl.table$IntronStart,ensembl.table$IntronEnd)
ensembl.table$intronID <- paste(ensembl.table$Chromosome, ensembl.table$cleanStart, ensembl.table$cleanEnd, sep = '_')

Granges.cassette.exons.left <- GRanges(seqnames = gsub(pattern = 'chr', replacement = '', cassette.exons$X.chrom),
                                  IRanges(start = cassette.exons$chromStart-50, end = cassette.exons$chromStart + 50))

Granges.cassette.exons.right <- GRanges(seqnames = gsub(pattern = 'chr', replacement = '', cassette.exons$X.chrom),
                                  IRanges(start = cassette.exons$chromEnd-50, end = cassette.exons$chromEnd + 50))



### Now look for overlap
introns.GRanges <- GRanges(seqnames = ensembl.table$Chromosome,
                           IRanges(start= ensembl.table$cleanStart, end = ensembl.table$cleanEnd))

overlap.right <- findOverlaps( Granges.cassette.exons.right, introns.GRanges)
overlap.left <- findOverlaps( Granges.cassette.exons.left, introns.GRanges)


#exons.with.long.introns.on.side <- intersect(overlap.right@queryHits, overlap.left@queryHits)

cassette.exons$intron.on.left <- overlap.left@subjectHits [ match(1:nrow(cassette.exons), table = overlap.left@queryHits) ]
cassette.exons$intron.on.right <- overlap.right@subjectHits [ match(1:nrow(cassette.exons), table = overlap.right@queryHits) ]

cassette.exons <-  subset(cassette.exons, ! is.na(intron.on.left) & !is.na(intron.on.right))  ### Now we take the ones with long introns on both side

cassette.exons$intronIDLeft <- ensembl.table$intronID[ cassette.exons$intron.on.left ]
cassette.exons$intronStartLeft <- ensembl.table$cleanStart[ cassette.exons$intron.on.left ]
cassette.exons$intronEndLeft <- ensembl.table$clceanEnd[ cassette.exons$intron.on.left ]
cassette.exons$intronLeftGene <- ensembl.table$Gene[ cassette.exons$intron.on.left ]
cassette.exons$intronLeftHUGO <- ensembl.table$HUGO[ cassette.exons$intron.on.left ]


cassette.exons$intronIDRight <- ensembl.table$intronID[ cassette.exons$intron.on.right ]
cassette.exons$intronStartRight <- ensembl.table$cleanStart[ cassette.exons$intron.on.right ]
cassette.exons$intronEndRight <- ensembl.table$cleanEnd[ cassette.exons$intron.on.right ]
cassette.exons$intronRightGene <- ensembl.table$Gene[ cassette.exons$intron.on.right ]
cassette.exons$intronRightHUGO <- ensembl.table$HUGO[ cassette.exons$intron.on.right ]




#### Now we have two cases, the first one being when the same intron overlaps both side of the cassette exon
cassette.exons.simple <-  subset(cassette.exons, intron.on.left == intron.on.right & cassette.exons$intronStartLeft < chromStart - minlength &  cassette.exons$intronEndRight > chromEnd + minlength)

#### sanity check with the NTM gene
print(subset(cassette.exons.simple, intronLeftHUGO ==  'NTM'))


cassette.exons.simple$junctionPos <- (cassette.exons.simple$chromStart +  cassette.exons.simple$chromEnd)/2
write.csv(x = cassette.exons.simple, file = 'data/cassette_exons_with_long_introns.csv', row.names = FALSE)
