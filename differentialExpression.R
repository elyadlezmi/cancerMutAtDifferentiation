### Libraries
library(edgeR)
library(fgsea)


### combine counts
test <- 'PSCs'
# test <- 'Differentiation_timeline'

gtf <- read.csv('/home/elyad/Bioinformatics/annotations/gencode.v38.Table.csv')[,c('ensemblID', 'gene_name')]
layout <- read.csv(paste(test,'layout.csv', sep="_"))

samples <- as.character(layout$Accession)
counts <- read.csv(paste0('/media/elyad/HDD/NeuralDiff/counts/', samples[1], 'ReadsPerGene.out.tab'), sep = '\t', skip = 1)[-c(1,2),-c(3,4)]
names(counts) <- c("ensemblID", samples[1]) 

for(sample in samples[2:length(samples)]) {
  counts[sample] <- read.csv(paste0('/media/elyad/HDD/NeuralDiff/counts/', sample, 'ReadsPerGene.out.tab'), sep = '\t', skip = 1)[-c(1,2),2] }

counts <- merge( gtf, counts, by='ensemblID')

write.csv(counts, paste(test, 'counts.csv', sep='_'))

### generate normalized count table 
counts <- read.csv(paste0(test,'_counts.csv'), row.names=1)
layout <- read.csv(paste0(test,'_layout.csv'))

y <- DGEList(counts[, -c(1,2)], genes=counts[c(1,2)])
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# calculate log2 transformed TMM normalized values
g <-  cpm(y, log=TRUE) 
g[g < -0.5] <- -0.5 # floor values to -0.5
write.csv(cbind(y$genes, g), paste0(test,'_normalized_expression.csv')) # write table

### DE analysis
layout$muts <- factor(layout$muts, levels = c("WT", "TP53"))
if(test == 'PSCs') {
    design <- model.matrix(~layout$CellLineID + layout$muts)}
if(test == 'Differentiation_timeline') {
    layout <- layout[layout$CellState == 'differentiated',] # remove PSCs
    design <- model.matrix(~layout$TimePoint + layout$muts)} 

counts <- counts[,c("ensemblID", "gene_name", layout$Accession)]
y <- DGEList(counts[, -c(1,2)], genes=counts[c(1,2)])
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
y <- glmFit(y, design)
y <- glmLRT(y)
DE <- data.frame(topTags(y,n = Inf,sort.by='p.value'))
write.csv(DE, paste0(test,"_DE.csv"))

### gsea analysis
df <- read.csv(paste0(test,'_DE.csv'), row.names = 1)
df$ranking <- -log(df$PValue) * sign(df$logFC)
write.table(df[order(df$ranking, decreasing = T), c("gene_name","ranking")], 'de.rnk', sep = '\t', row.names = F, col.names = F, quote = F)

ranks <- read.table('de.rnk',header=FALSE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$V2,ranks$V1)

gene_sets <- c(gmtPathways('c2.cp.kegg.v7.5.1.symbols.gmt'),gmtPathways('c5.go.bp.v7.5.1.symbols.gmt'))

fgseaRes <- fgsea(gene_sets, ranks, eps=0, minSize=50)
fgseaRes <- fgseaRes[padj <= 0.05][order(pval)]
collapsed <- collapsePathways(fgseaRes, gene_sets, ranks)
fgseaRes <- fgseaRes[pathway %in% collapsed$mainPathways][order(pval)]

fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)

write.csv(fgseaRes, paste0(test,'_gsea.csv'))

