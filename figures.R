
################## Libraries #################

library(gtools) # install.packages("OpenRepGrid")
library(ggplot2)
library(pals)
library(OpenRepGrid)
library(plyr)
library(dplyr)
library(Polychrome)
library(ggpattern) # install.packages("ggpattern")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr)
library(pathview)
pal.bands(stepped()) # build-in color palette

DIR <- './Results/Tables/'

######### Fig 1 --------------------------------------------
### A ----
# ES/iPS - inner circle
df <- read.csv(file = paste0(DIR,'finalMetaData.csv'))

col <- df$iPSES
p <- data.frame(table(col))
p <- p[order(p$Freq, decreasing = T),]
p$col <- factor(p$col, levels=p$col)
p$perc <- p$Freq / sum(p$Freq)
print(head(p, n = 10))
print(length(unique(col)))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect(color="black") + scale_fill_manual(values=recycle(rev(stepped())[c(22, 18)], length(p$col), na.fill = FALSE)) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none")

# Cell-lines - outer circle
col <- df[df$iPSES == 'ES',]$CellLineID
p1 <- data.frame(table(col))
p1 <- p1[order(p1$Freq, decreasing = T),]

col <- df[df$iPSES == 'iPS',]$CellLineID
p2 <- data.frame(table(col))
p2 <- p2[order(p2$Freq, decreasing = T),]

p <- rbind(p1, p2)
p$perc <- p$Freq / sum(p$Freq)
p$col[duplicated(p$col)]
p$col <- factor(p$col, levels=p$col)
print(head(p, n = 10))
print(length(unique(p$col)))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect(color="black") + scale_fill_manual(values=recycle(rev(stepped())[c(9:24, 1:8)], length(p$col), na.fill = FALSE)) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none") 

### B ----
# diff type by germlayer - inner circle
df[df$CellState == 'pluripotent', 'Lineage'] <- 'undifferentiated'
df[df$CellState == 'pluripotent', 'DiffType'] <- 'undifferentiated'

col <- df[df$Lineage == 'Ectoderm',]$DiffType
p1 <- data.frame(table(col))
p1 <- p1[order(p1$Freq, decreasing = T),]

col <- df[df$Lineage == 'Mesoderm',]$DiffType
p2 <- data.frame(table(col))
p2 <- p2[order(p2$Freq, decreasing = T),]

col <- df[df$Lineage == 'Endoderm',]$DiffType
p3 <- data.frame(table(col))
p3 <- p3[order(p3$Freq, decreasing = T),]

col <- df[df$Lineage == 'undifferentiated',]$DiffType
p4 <- data.frame(table(col))
p4 <- p4[order(p4$Freq, decreasing = T),]

p <- rbind(p1, p2, p3, p4)
p$perc <- p$Freq / sum(p$Freq)
p$col <- factor(p$col, levels=p$col)
p$col[duplicated(p$col)]
print(head(p, n = 10))
print(length(unique(col)))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect(color="black") + scale_fill_manual(values=recycle(stepped(), length(p$col), na.fill = FALSE)) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none") 

# germ layer - outer circle
col <- df$Lineage
p <- data.frame(table(col))[c(1,3,2,4),]
p$perc <- p$Freq / sum(p$Freq)
p$col <- factor(p$col, levels=c('Ectoderm','Mesoderm', 'Endoderm', 'undifferentiated'))
print(head(p, n = 10))
print(length(unique(col)))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect(color="black") + scale_fill_manual(values=stepped()[c(15, 3, 7, 12)]) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none")

### C - right side ----
# pie muts in cell-lines - upper
df <- read.csv(file = paste0(DIR,'finalMetaData.csv'))
df$is_mut <- df$n_mut > 0

a <- length(unique(df$CellLineID))
df$mutCellLineID <- paste0(df$CellLineID, as.character(df$is_mut))
b <- length(unique(df$mutCellLineID)) - a
p <- data.frame(col1=c('Mutated Cell-Lines', 'Clean Cell-Lines'),col2=c(b, a-b))

ggplot(p, aes(x="", y=col2)) +
  geom_bar_pattern(stat='identity', width=1, pattern = c('stripe', 'none'), fill=stepped()[18],
                   color='black', pattern_fill=stepped()[21], pattern_density=0.25,
                   pattern_spacing=0.04, pattern_angle=45) +
  coord_polar("y", start=0) + theme_void() + theme(legend.position="none")
                       
# pie muts in studies - lower
a <- length(unique(df$Study))
df$mutStudy <- paste0(df$Study, as.character(df$is_mut))
b <- length(unique(df$mutStudy)) - a
p <- data.frame(col1=c('Mutated Studies', 'Clean Studies'),col2=c(b, a-b))

ggplot(data=p, aes(x="", y=col2)) +
  geom_bar_pattern(stat='identity', width=1, pattern = c('stripe', 'none'), fill=stepped()[19],
                   color='black', pattern_fill=stepped()[21], pattern_density=0.25,
                   pattern_spacing=0.04, pattern_angle=45) +
  coord_polar("y", start=0) + theme_void() + theme(legend.position="none")

##### for Figs C-E ----
mutations <- read.csv(paste0(DIR,'Mutations.csv'), skip = 1)
mutations$acquisition <- "Unknown (possible germline)"
mutations[mutations$is_somatic == "True", "acquisition"] <- "Acquired in vitro (undifferentiated)"
recentlyAcquired <- read.csv(paste0(DIR,'AcquiredMutations.csv'))
recentlyAcquired <- recentlyAcquired$Accession 
mutations[mutations$Accession %in% recentlyAcquired, "acquisition"] <- "Acquired / expanded during differentiation"
mutations$acquisition <- 
  factor(mutations$acquisition, levels=c("Acquired in vitro (undifferentiated)", "Unknown (possible germline)", 
                                         "Acquired / expanded during differentiation"))

df <- read.csv(file = paste0(DIR,'finalMetaData.csv'))
df[df$CellState == 'pluripotent', 'Lineage'] <- 'Undifferentiated'
df$is_mut <- df$n_mut > 0 
# keep only genes mutated in >= 2 samples
mutatedIn2 <- read.csv(paste0(DIR,'MutationSummarizedPerGene.csv'))
mutatedIn2 <- mutatedIn2[mutatedIn2$N_samples > 1, 'Gene']
mutations <- mutations[mutations$Gene %in% mutatedIn2, ]

mut_order <- table(mutations$Gene)
mut_order <- mut_order[order(mut_order, decreasing = T)]
mutations$Gene <- factor(mutations$Gene, levels = names(mut_order))

### C - left side ----
# circle of muts to all by ES/iPS - upper
df$newiPSES <- paste0(df$iPSES, as.character(df$is_mut))
col <- df$newiPSES
p <- data.frame(table(col))[c(1, 2, 3, 4),]
p$perc <- p$Freq / sum(p$Freq)
p$col <- factor(p$col, levels=c('ESFALSE', 'ESTRUE', 'iPSFALSE', 'iPSTRUE'))
print(head(p, n = 10))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect_pattern(color='black', pattern = rep(c('none', 'stripe'),2), 
                    pattern_fill=stepped()[21], pattern_density=0.5, pattern_spacing=0.02, pattern_angle=45) +
  scale_fill_manual(values=stepped()[c(3, 3, 7, 7)]) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none")

# circle of muts to all by lineage - lower
df$newLineage <- paste0(df$Lineage, as.character(df$is_mut))
col <- df$newLineage
p <- data.frame(table(col))[c(1, 2, 5, 6, 3, 4, 7, 8),]
p$perc <- p$Freq / sum(p$Freq)
p$col <- factor(p$col, levels=c('EctodermFALSE', 'EctodermTRUE', 'MesodermFALSE', 'MesodermTRUE', 
                                'EndodermFALSE', 'EndodermTRUE', 'UndifferentiatedFALSE', 'UndifferentiatedTRUE'))
print(head(p, n = 10))
p$ymax = cumsum(p$perc)
p$ymin = c(0, head(p$ymax, n=-1))

ggplot(data=p, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=col)) +
  geom_rect_pattern(color='black', pattern = rep(c('none', 'stripe'),4), 
                    pattern_fill=stepped()[21], pattern_density=0.5, pattern_spacing=0.02, pattern_angle=45) +
  scale_fill_manual(values=stepped()[c(15, 15, 3, 3, 7, 7, 12, 12)]) +
  coord_polar(theta="y") + xlim(c(2, 4)) + theme_void() + theme(legend.position="none")

### D top ----
# mutated genes bars by CellState - left side
p <- data.frame(table(mutations$Gene, mutations$CellState))
names(p) <- c("Gene", "var", "Samples")
totalSamples <- aggregate(p$Samples, by=list(p$var), FUN=sum)
totalSamples <- setNames(totalSamples[,2], totalSamples[,1])   
for(i in names(totalSamples)){
  p[p$var == i, "Samples"] <- (p[p$var == i, "Samples"]/totalSamples[i])}

ggplot(data=p, aes_string(x="var", y="Samples", fill="Gene")) + ylab("Proprtion of mutated samples") +
  geom_bar(color='black', stat="identity") +
  scale_fill_manual(values=c(stepped()[1:20], "grey90", "grey75", "grey60", "grey45", "grey30", "grey15", "grey5")) +
  theme_classic() + theme(legend.position="none", axis.title.x=element_blank()) + scale_y_continuous(expand = c(0,0))

# mutated genes bars by cell origin - right side
p <- data.frame(table(mutations$Gene, mutations$iPSES))
names(p) <- c("Gene", "var", "Samples")
totalSamples <- aggregate(p$Samples, by=list(p$var), FUN=sum)
totalSamples <- setNames(totalSamples[,2], totalSamples[,1])   
for(i in names(totalSamples)){
  p[p$var == i, "Samples"] <- (p[p$var == i, "Samples"]/totalSamples[i])}

ggplot(data=p, aes_string(x="var", y="Samples", fill="Gene")) + ylab("Proprtion of mutated samples") +
  geom_bar(color='black', stat="identity") +
  scale_fill_manual(values=c(stepped()[1:20], "grey90", "grey75", "grey60", "grey45", "grey30", "grey15", "grey5")) +
  theme_classic() + theme(legend.position="none", axis.title.x=element_blank()) + scale_y_continuous(expand = c(0,0))

### D bottom ----
# ES\iPS acquired in vitro - left side
p <- mutations[mutations$acquisition != "Unknown (possible germline)",]
p <- data.frame(table(p$Gene, p$iPSES))
names(p) <- c("Gene", "iPSES", "Samples")

ggplot(data=p, aes_string(x="iPSES", y="Samples", fill="Gene")) + 
  ylab("Number of mutations acquired in vitro") +
  geom_bar(color='black', stat="identity") + 
  geom_hline(yintercept = 25, color=stepped()[23]) +
  scale_fill_manual(values=c(stepped()[1:20], "grey90", "grey75", "grey60", "grey45", "grey30", "grey15", "grey5")) +
  theme_classic() + theme(legend.position="none", axis.title.x=element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,330))

# ES\iPS germline\unknown - right side
p <- mutations[mutations$acquisition == "Unknown (possible germline)",]
p <- data.frame(table(p$Gene, p$iPSES))
names(p) <- c("Gene", "iPSES", "Samples")

ggplot(data=p, aes(x=iPSES, y=Samples, fill=Gene)) + 
  ylab("Number of mutations of unknown origin") +
  geom_bar(color='black', stat="identity") + 
  geom_hline(yintercept = 25, color=stepped()[23]) + 
  scale_fill_manual(values=c(stepped()[1:20], "grey90", "grey75", "grey60", "grey45", "grey30", "grey15", "grey5")) +
  theme_classic() + theme(legend.position="none", axis.title.x=element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,330))

### E ----
# AF density per mut origin - upper
ggplot() +
  geom_density(data=mutations[mutations$acquisition == 'Acquired in vitro (undifferentiated)',], aes(x=AF,y=0.1*..count..),
               size=1.5, bw=0.1, alpha=0.2, color=stepped()[15], fill=stepped()[15]) +
  geom_density(data=mutations[mutations$acquisition == 'Unknown (possible germline)',], aes(x=AF,y=0.1*..count..),
               size=1.5, bw=0.1, alpha=0.2, color=stepped()[7], fill=stepped()[7]) +
  geom_density(data=mutations[mutations$acquisition == 'Acquired / expanded during differentiation',], aes(x=AF,y=0.1*..count..),
               size=1.5, bw=0.1, alpha=0.2, color=stepped()[3], fill=stepped()[3]) +
  theme_classic() + xlab('Allelic Frequency') + ylab('Number of mutations') + 
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0), limits=(c(0,1))) + 
  geom_vline(xintercept=0.5, linetype="dotted")

# AF density by TP53 to others - lower
genes <- read.csv(paste0(DIR,'MutationSummarizedPerGene.csv'))
genes <- genes[genes$N_lines >= 2, "Gene"]

ggplot(data=mutations[mutations$Gene %in% genes,], aes(x=AF,y=0.1*..count..)) +
  geom_density(aes(color=Gene=='TP53', fill=Gene=='TP53'), size=1.5, bw=0.1, alpha=0.2) +
  theme_classic() + xlab('Allelic Frequency') + ylab('Number of mutations') +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0), limits=(c(0,1)))  + 
  geom_vline(xintercept=0.5, linetype="dotted") + 
  scale_color_manual(values=c(stepped()[c(22,1)])) +
  scale_fill_manual(values=c(stepped()[c(22,1)]))

### F&G ----
# mapped P53 mutations appearances - fig F
df <- read.csv(file = paste0(DIR, 'Mutations.csv'), skip = 1)

p <- df[df$Gene == "TP53",]
p$AA_POS <- as.numeric(sub("\\D*(\\d+).*", "\\1", p$Mutation.AA))
p <- p  %>% group_by(COSMIC_ID) %>% dplyr::summarise(
  Accession = length(Accession), POS=first(AA_POS), Mutation.AA=first(Mutation.AA))
p$Mutation.AA <- substring(p$Mutation.AA, 3)

ggplot(p, aes(x=POS, y=Accession, label = Mutation.AA)) +
  geom_segment( aes(x=POS, xend=POS, y=0, yend=Accession), size=1.3, alpha=0.9, color="black") +
  geom_label_repel(aes(label=Mutation.AA), box.padding = 0.5, max.overlaps = Inf) +
  theme_classic() + xlim(c(0,395)) + xlab("Position along p53 (a.a)") +
  ylab("Number of Samples per Mutation") + scale_y_continuous(expand = c(0,0), limits = c(0,130))

# comparison to cosmic tumors - fig G
t <- p
df <- read.csv(file = paste0(DIR, 'V96_38_MUTANTCENSUS.csv'))

p <- df[!df$GENOMIC_MUTATION_ID %in% c("null"),]
p$POS <- as.numeric(sub("\\D*(\\d+).*", "\\1", p$MUTATION_AA))
p <- p  %>% group_by(GENOMIC_MUTATION_ID) %>% dplyr::summarise(
  POS=first(POS), num=length(GENOMIC_MUTATION_ID), COSMIC_ID=first(GENOMIC_MUTATION_ID))
p$WE_HAVE <- p$COSMIC_ID %in% t$COSMIC_ID
p <- p  %>% group_by(POS) %>% dplyr::summarise(
  num=sum(num), WE_HAVE=sum(WE_HAVE))
p$WE_HAVE_BOOL <- as.character(p$WE_HAVE > 0)

ggplot(p, aes(x=POS, y=num)) +
  geom_segment(aes(x=POS, xend=POS, y=0, yend=num, color=WE_HAVE_BOOL), size=1.3, alpha=0.9) +
  theme_classic() +
  theme(legend.position = "none", panel.border = element_blank(),) + 
  xlab("Position along p53 (a.a)") +
  ylab("Number of Substitutions in COSMIC (x1000)") + scale_color_manual(values=c('TRUE' = "BLACK", 'FALSE' = "DARKGRAY" )) +
  scale_y_continuous(expand = c(0,0))

######### Fig 2 --------------------------------------------

### A ----
#  plot gsea (ESCs)
df <- read.csv('PSCs_gsea.csv')
df <- df[df$NES < 0,]
df <- head(df,n=5)
df$pathway <- factor(df$pathway, levels=rev(df$pathway))

ggplot(data=df, aes(x=pathway, y=-log(padj))) +
    geom_bar(stat = "identity", fill=stepped()[13], color="black", alpha=0.7) + 
    geom_hline(yintercept=-log(0.05), color=stepped()[2]) +
    theme_classic() + theme_classic() + coord_flip() +
    scale_y_continuous(expand = c(0,0))

### B ----
# plot kegg 
df <- read.csv('PSCs_DE.csv', row.names = 1)[c("gene_name","logFC","FDR")]
df <- df[order(abs(df$logFC), decreasing = T),]
df <- df[!duplicated(df$gene_name),]
row.names(df) <- df[,1]
df <- data.matrix(df)
df[,2] <- as.numeric(df[,2])
pathview(gene.data = df[,2], pathway.id = "hsa04115",species = "hsa", out.suffix = 'PSCs', gene.idtype="SYMBOL", limit=list(gene=2), 
         low=list(gene =  stepped()[13], cpd = "dodgerblue4"), high=list(gene = stepped()[1], cpd = "red2"), kegg.native=F)

### C ----
#  plot expression of MDM2 in ESCs 
gene <- 'MDM2'

df <- read.csv('PSCs_normalized_expression.csv', row.names = 1)
df <- df[df$gene_name == gene, -c(1,2)] # CDKN1A, MDM2

df <- reshape2::melt(df)
layout <- read.csv('PSCs_layout.csv')

df$mut <- layout$muts
df$mut <- factor(df$mut, levels=c("WT", "TP53"))
df$CellLineID <- layout$CellLineID

for(cl in unique(df$CellLineID)) {
    
    for(c in c('TP53', 'WT')) {
        samples <- df[df$CellLineID == cl & df$mut == c, 'variable']
        values <- df[df$CellLineID == cl & df$mut == c, 'value']
        
        top <- quantile(values, 0.75) + 1.5*IQR(values); bottom <- quantile(values, 0.25) - 1.5*IQR(values)
        outliers <- samples[values <= bottom | values >= top]
        
        df <- df[!df$variable %in% outliers,] }
    
    p <- df[df$CellLineID == cl,]
    print(cl)
    print(wilcox.test(p[p$mut == 'WT', "value"], p[p$mut == 'TP53', "value"]))  }


ggplot(data=df, aes(x=mut, y=value, fill=mut)) + 
    geom_boxplot(color="black", coef=50, width=0.5) + 
    geom_jitter(width=0.1, color="black", alpha=0.5) +
    # geom_bar(color="black", stat="summary", position=position_dodge(), width=0.5, alpha=1) +
    # stat_summary(fun.data = mean_se, geom = "errorbar", size=0.25, width=.1,position=position_dodge(), color="black") +
    theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    facet_wrap(vars(CellLineID)) + ylab(paste(gene, 'expression', sep= ' ')) + scale_fill_manual(values=stepped()[c(6,2)]) + 
    scale_y_continuous(expand = c(0,0)) 

### D ----
#  plot gsea (differentiation)
df <- read.csv('Differentiation_timeline_gsea.csv')
df <- df[df$NES < 0,]
df <- head(df,n=5)
df$pathway <- factor(df$pathway, levels=rev(df$pathway))

ggplot(data=df, aes(x=pathway, y=-log(padj))) +
    geom_bar(stat = "identity", fill=stepped()[13], color="black", alpha=0.7) + 
    geom_hline(yintercept=-log(0.05), color=stepped()[2]) +
    theme_classic() + theme_classic() + coord_flip() +
    scale_y_continuous(expand = c(0,0))

### E & F ----
tmm <- read.csv('Differentiation_timeline_normalized_expression.csv', row.names = 1)
layout <- read.csv('Differentiation_timeline_layout.csv')

gene <-  'POU5F1'
# gene <-  'MAP2'

df <- tmm[tmm$gene_name == gene, -c(1,2)] # CDKN1A, MDM2
df <- reshape2::melt(df)
df$mut <- layout$muts
df$mut <- factor(df$mut, levels=c('WT', 'TP53'))
df$TimePoint <- layout$TimePoint
df <- df[df$TimePoint %in% c(0,5,7),]

for(tp in unique(df$TimePoint)) {
    for(c in c('TP53', 'WT')) {
        samples <- df[df$TimePoint == tp & df$mut == c, 'variable']
        values <- df[df$TimePoint == tp & df$mut == c, 'value']
        
        top <- quantile(values, 0.75) + 1.5*IQR(values); bottom <- quantile(values, 0.25) - 1.5*IQR(values)
        outliers <- samples[values < bottom | values > top]
        
        df <- df[!df$variable %in% outliers,] } 
    
    p <- df[df$TimePoint == tp,]
    print(tp)
    print(wilcox.test(p[p$mut == 'WT', "value"], p[p$mut == 'TP53', "value"])) }

ggplot() +
    geom_smooth(data=df[df$mut=='WT',],aes(x=TimePoint, y=value), color=stepped()[6], se = F, method = "lm") + 
    geom_smooth(data=df[df$mut!='WT',],aes(x=TimePoint, y=value), color=stepped()[2], se = F, method = "lm") + 
    geom_boxplot(data=df[df$mut=='WT',],aes(x=TimePoint-0.15, y=value, group=TimePoint), coef=50, width=0.3, fill=stepped()[6]) + 
    geom_boxplot(data=df[df$mut!='WT',],aes(x=TimePoint+0.15, y=value, group=TimePoint), coef=50, width=0.3, fill=stepped()[2]) + 
    geom_point(data=df[df$mut=='WT',],aes(x=TimePoint-0.15, y=value), alpha=0.2) + 
    geom_point(data=df[df$mut!='WT',],aes(x=TimePoint+0.15, y=value), alpha=0.2) + 
    scale_fill_manual(values=stepped()[c(6,2)]) + theme_classic() + ylab(paste(gene, 'expression', sep=' ')) +
    scale_y_continuous(expand = c(0,0)) + xlab('Days into neural differentiation')

