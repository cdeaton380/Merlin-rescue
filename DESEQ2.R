library("gdata")
library(DESeq2)
library(writexl)
library(ggplot2)
library(org.Hs.eg.db)




########################## Count Matrix Input ################################

pasCts <- read.delim("featurecounts_pc.Rmatrix.txt",sep="\t",
                     header = T,stringsAsFactors = F,row.names = 1)

cts <- as.matrix(pasCts)



# Metadata
MetaData=data.frame(Type=colnames(cts
))



MetaData$Group=c(rep('dox',times=4),rep('no_dox',times=3))


rownames(MetaData)=colnames(cts)

colnames(MetaData)[2]="Group"

MetaData$Type<- factor(MetaData$Type)
MetaData$Group<- factor(MetaData$Group)


########################## Differential Expression Analisys ################################

# Contruct DeseqDataSet

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = MetaData,
                              design = ~ Group)

rownames(dds)
colnames(dds)

dds <- DESeq(dds)

colData(dds)


# Filtering by count

nrow(dds)  
dds_filt=dds[rowSums(counts(dds)) >=10,]

nrow(dds_filt)  

dds=dds_filt  

dds <- DESeq(dds)

nrow(dds) 

# Storing the filtered and normalized DESeq2 Object

saveRDS(dds,"normalizedDESeq2object.rds")



# Transforming count data 

vsd <- vst(dds, blind=FALSE)

vsd_transformed_data=assay(vsd)
vsd_transformed_data=as.data.frame(vsd_transformed_data)
vsd_transformed_data$Genes=rownames(vsd_transformed_data)

vsd_transformed_data=vsd_transformed_data[,c(8,1:7)]

writexl::write_xlsx(vsd_transformed_data,"vst_transformed.xlsx")



# Heatmap based on top 2000 most variable genes

topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),2000)
vsd_data_heatmap=vsd_transformed_data[topVarGenes,]

vsd_data_heatmap_matrix=as.matrix(vsd_data_heatmap[,2:8])


library(pheatmap)

jpeg("Heatmap_SCALED.jpeg",width=11,height = 8,res=190,units="in")
pheatmap(vsd_data_heatmap_matrix, annotation_col = MetaData,scale='row',
         show_rownames = F, show_colnames = T)

dev.off()


# PCA Plot

jpeg("PCA.jpeg",width=9,height=8,res=150,units="in")
p=plotPCA(vsd,intgroup="Type")

p

dev.off()



################### Fold Change 
res <- 
  results(dds, contrast=c("Group","dox","no_dox"),alpha=0.05) # CHANGE 

head(res)

res=as.data.frame(res)

res$Genes=rownames(res)

res=res[,c(7,1:6)]

summary(res)


i=which(is.na(res$symbol))
res$symbol[i]=res$Genes[i]

write_xlsx(res,"All_DEG_dox_vs_no_dox.xlsx") 



# Volcano Plot 

alpha <- 0.05 # Threshold on the p-value

jpeg("Volcano_Plot_dox_vs_no_dox.jpeg",width=15,height=12,
     res=150,units="in") 

with(res,plot(log2FoldChange, -log10(res$pvalue), 
              col="grey", panel.first=grid(),
              main="dox vs no dox", xlab="log2(fold-change)", 
              ylab="-log10(p-value)",
              pch=20, cex=0.8, xlim=c(-11,10)) )  
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

with(subset(res, pvalue<0.05 & 
              log2FoldChange>=1), points(log2FoldChange, -log10(pvalue), 
                                         pch=20, col="red", cex=0.8))
with(subset(res, pvalue<0.05 & 
              log2FoldChange<=-1), points(log2FoldChange, -log10(pvalue), 
                                          pch=20, col="blue", cex=0.8))

legend("topright",legend=c("Upregulated","Downregulated","Non-Significant"),
       pch=20,col=c("red","blue","grey"),cex=0.9)

genes_selected=res[order(res$pvalue),]
genes_selected=genes_selected[1:10,] # Getting top 10 highly significant genes based on pValue

sig=subset(res,pvalue<0.05)
sig=sig[order(sig$log2FoldChange,decreasing = T),]

genes_selected=rbind(genes_selected,sig[1:5,],tail(sig,n=5))
genes_selected=genes_selected[!duplicated(genes_selected$Genes),]


library(basicPlotteR)
addTextLabels(genes_selected$log2FoldChange,
              -log10(genes_selected$pvalue), genes_selected$Genes, 
              cex.label=1, col.label="black")

dev.off()









