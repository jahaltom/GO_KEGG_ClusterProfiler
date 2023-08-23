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
BiocManager::install(organism, character.only = TRUE,ask=FALSE)
library(organism, character.only = TRUE)


fileName <- "OrigList"
conn <- file(fileName,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
    f=linn[i]




    # reading in data from deseq2
    df = read.csv(paste(f,"/Differential_expression_analysis_table.csv",sep=""), header=TRUE)
    #Subset sig
    #df=df[ which(df$pvalue <= 0.05), ]
    
    #Add gene names. Will be longer than th eone from org.Hs.eg.db. 
    gn = read.csv("GeneAnn.txt",sep='\t',header=FALSE)
    colnames(gn)=c("ID","SYMBOL")
    dfAnn=unique((merge(gn,df,by="ID",all.y=TRUE)))
    write.table(dfAnn,paste("out/",f,"/Differential_expression_analysis_table.GeneName.csv",sep=""),sep = '\t',row.names = FALSE)
    

    
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
    
    # Produce the 1st  native KEGG plot (PNG)
     
    dme <-pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[1], species = kegg_organism)
 
    
    # Produce the 2nd  native KEGG plot (PNG)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[2], species = kegg_organism)
    
    #A bunch of mis KEGG plots
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa01100", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa01230",out.suffix="Biosynthesis of amino acids", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa01232",out.suffix="Nucleotide metabolism", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00010",out.suffix="Glycolysis Gluconeogenesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00020",out.suffix="Citrate cycle (TCA cycle)", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00520",out.suffix="Amino sugar and nucleotide sugar metabolism", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00620",out.suffix="Pyruvate metabolism", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00061",out.suffix="Fatty acid biosynthesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00062",out.suffix="Fatty acid elongation", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00071",out.suffix="Fatty acid degradation", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa01040",out.suffix="Biosynthesis of unsaturated fatty acids", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00230",out.suffix="Purine metabolism", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00240",out.suffix="Pyrimidine metabolism", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00290",out.suffix="Valine, leucine and isoleucine biosynthesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00300",out.suffix="Lysine biosynthesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00400",out.suffix="Phenylalanine, tyrosine and tryptophan biosynthesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00510",out.suffix="N-Glycan biosynthesis", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00513",out.suffix="Various types of N-glycan biosynthesis", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00512",out.suffix="Mucin type O-glycan biosynthesis", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00515",out.suffix="Mannose type O-glycan biosynthesi", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03020",out.suffix="RNA polymerase", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03022",out.suffix="Basal transcription factors", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03040",out.suffix="Spliceosome", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03010",out.suffix="Ribosome", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00970",out.suffix="Aminoacyl-tRNA biosynthesis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03013",out.suffix="Nucleocytoplasmic transport", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03015",out.suffix="mRNA surveillance pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03008",out.suffix="Ribosome biogenesis in eukaryotes", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03060",out.suffix="Protein export", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04141",out.suffix="Protein processing in endoplasmic reticulum", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03262",out.suffix="Virion - Coronavirus New!", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa02010",out.suffix="ABC transporters", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa02060",out.suffix="Phosphotransferase system (PTS)", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04010",out.suffix="MAPK signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04630",out.suffix="JAK-STAT signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04064",out.suffix="NF-kappa B signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04668",out.suffix="TNF signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04066",out.suffix="HIF-1 signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04151",out.suffix="PI3K-Akt signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04152",out.suffix="AMPK signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04150",out.suffix="mTOR signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04061",out.suffix="Viral protein interaction with cytokine and cytokine receptor", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04144",out.suffix="Endocytosis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04145",out.suffix="Phagosome", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04142",out.suffix="Lysosome", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04146",out.suffix="Peroxisome", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04140",out.suffix="Autophagy - animal", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04137",out.suffix="Mitophagy - animal", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04210",out.suffix="Apoptosis", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04620",out.suffix="Toll-like receptor signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04621",out.suffix="NOD-like receptor signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04622",out.suffix="RIG-I-like receptor signaling pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04623",out.suffix="Cytosolic DNA-sensing pathway", species = kegg_organism)
    dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05171",out.suffix="Coronavirus disease - COVID-19", species = kegg_organism)


    
    
    
    
    
    files=list.files(pattern = "*hsa*")
    for (i in files){
        file.move(i, paste("out/",f,"/",sep=""))

    }

}
close(conn)   
