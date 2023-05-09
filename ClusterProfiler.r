library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggridges)
library(ggupset)
library(pathview)
library(filesstrings)
########GO#########

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


fileName <- "list"
conn <- file(fileName,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
    f=linn[i]




    # reading in data from deseq2
    df = read.csv(paste(f,"/Differential_expression_analysis_table.csv",sep=""), header=TRUE)
    #Subset sig
    #df=df[ which(df$pvalue <= 0.05), ]
    
    
    
    # we want the log2 fold change 
    original_gene_list <- df$log2FoldChange
    
    # name the vector gene id
    names(original_gene_list) <- df$ID
    

    
    # Convert gene IDs for gseKEGG function
    # We will lose some genes here because not all IDs will be converted
    ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=organism)
     # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
    
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    df2 = df[df$ID %in% dedup_ids$ENSEMBL,]
    
    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$GeneName = dedup_ids$SYMBOL
    
    # Create a vector of the gene unuiverse
    gene_list <- df2$log2FoldChange
    
    # Name vector with ENTREZ ids
    names(gene_list) <- df2$GeneName
    
    # omit any NA values 
    gene_list<-na.omit(gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)
    
    write.table(df2,paste("out/",f,"/Differential_expression_analysis_table.GeneName.csv",sep=""),sep = '\t',row.names = FALSE)

    
    
    
    
    
    
    GOFunc <- function(GO)    
    {    
     
        gse <- gseGO(geneList=gene_list, 
                     ont =GO, 
                     keyType = "SYMBOL", 
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = organism, 
                     pAdjustMethod = "none")
        
        require(DOSE)
        
        #Dotplot
        pdf(paste("out/",f,"/GO_",GO,".Dotplot.pdf",sep=""),width=10,height=13)    
        print(dotplot(gse, showCategory=10, split=".sign",font.size = 11) + facet_grid(.~.sign))    
        dev.off()    
        
        #Encrichment Map
        pdf(paste("out/",f,"/GO_",GO,".EncrichmentMap.pdf",sep=""),width=10,height=13)    
        print(emapplot(pairwise_termsim(gse), showCategory = 10,cex_label_category = 0.85))    
        dev.off()    
        
        
        #Category Netplot
        # categorySize can be either 'pvalue' or 'geneNum'
        pdf(paste("out/",f,"/GO_",GO,".Netplot.pdf",sep=""),width=20,height=25)    
        print(cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3,font.size = 7))    
        dev.off()    
        
        #Ridgeplot
        pdf(paste("out/",f,"/GO_",GO,".Ridgeplot.pdf",sep=""),width=15,height=20)    
        print(ridgeplot(gse) + labs(x = "enrichment distribution",font.size = 7))    
        dev.off()       
        
        
        
      
        
        
    }  
    GOFunc("BP")
    GOFunc("CC")
    GOFunc("MF")
    
    
    
    ##########KEGG##########
    
    
    # Convert gene IDs for gseKEGG function
    # We will lose some genes here because not all IDs will be converted
    ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
     # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
    
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    df2 = df[df$ID %in% dedup_ids$ENSEMBL,]
    
    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$Y = dedup_ids$ENTREZID
    
    # Create a vector of the gene unuiverse
    kegg_gene_list <- df2$log2FoldChange
    
    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df2$Y
    
    # omit any NA values 
    kegg_gene_list<-na.omit(kegg_gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
              
    kegg_organism = "hsa"
    kk2 <- gseKEGG(geneList     = kegg_gene_list,
                   organism     = kegg_organism,
                   nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "none",
                   keyType       = "ncbi-geneid")
    
    
    #Dotplot
    pdf(paste("out/",f,"/KEGG_Dotplot.pdf",sep=""),width=10,height=13)    
    print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign",font.size = 11) + facet_grid(.~.sign))    
    dev.off()  
    
    
    #Encrichment Map
    pdf(paste("out/",f,"/KEGG_EncrichmentMap.pdf",sep=""),width=10,height=13)  
    print(emapplot(pairwise_termsim(kk2), showCategory = 10,cex_label_category = 0.85))    
    dev.off()   
    
    #Ridgeplot
    pdf(paste("out/",f,"/KEGG_Ridgeplot.pdf",sep=""),width=15,height=20)    
    print(ridgeplot(kk2) + labs(x = "enrichment distribution",font.size = 7))    
    dev.off()    
    
      
    #All KEGG Pathways
    #https://www.genome.jp/kegg/pathway.html
    
    # Produce the native KEGG plot (PNG)
     
    dme <-pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[1], species = kegg_organism)
 
    
    # Produce the native KEGG plot (PNG)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[2], species = kegg_organism)
    
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00190", species = kegg_organism)
    
    files=list.files(pattern = "*hsa*")
    for (i in files){
        file.move(i, paste("out/",f,"/",sep=""))

    }

}
close(conn)   







