library(clusterProfiler)
library(DESeq2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)



counts = read.table("raw_counts.csv",header=TRUE,sep = ',',check.names=FALSE,row.names=1)
#Add gene ID column to counts
counts=cbind(ID = rownames(counts), counts)
 # vector gene id
original_gene_list <- rownames(counts)


#Get gene name from Gen ID. Not all will have gene names, these are dropped.
ids<-bitr(original_gene_list, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=organism)
# remove duplicate IDS keeping 1st
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]


#merge
counts = counts[counts$ID %in% dedup_ids$ENSEMBL,]

# Create a new column  with the corresponding gene name
counts$GeneName = dedup_ids$SYMBOL


write.table(counts,paste("Differential_expression_analysis_table.AKvsZ.csv",sep=""),sep = '\t',row.names = FALSE)





#Read in metadata
md = read.table("metadata.tsv",header=TRUE,sep = '\t')
rownames(md)=md$SampleID
#Reorder counts cols to match md
counts=counts[,rownames(md)]
#Make DeSeq object
dds <- DESeqDataSetFromMatrix(counts,colData= md, design = ~Batch + Condition)



##DESeq2
dds = DESeq(dds)

#Output normalized counts
normCounts=counts(dds, normalized=TRUE)
normCounts = cbind(ID = rownames(normCounts), normCounts)
write.table(normCounts,"NormalizedCounts.tsv" ,sep = '\t',row.names = FALSE)

vsd <- vst(dds, blind=FALSE)
vsd = assay(vsd)
vsd = cbind(ID = rownames(vsd), vsd)
write.table(vsd,"VSDCounts.tsv" ,sep = '\t',row.names = FALSE)




list = read.table("Combo.tsv",sep = '\t',quote="")

for (x in 1:nrow(list)) {
    #########
    #Contrast
    result = results(dds, contrast=c("Condition",list[x,1],list[x,2]))
    ## Remove rows with NA
    result = result[complete.cases(result),]
    #Put GeneID as column
    result = cbind(ID = rownames(result), result)
    write.table(result,paste(list[x,1],"vs",list[x,2],".DGE.tsv",sep='_') ,sep = '\t',row.names = FALSE)

}
