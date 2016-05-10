####################################### As the Genes Do not contain information on survival times, we look for pathways ###################
########################################## As the PC doesnot contain much survival we use cox scores #######################################

library(globaltest)
library(org.Hs.eg.db)
library(GO.db)
library(KEGG.db)


exprs.norm <- ExpressionSet(Y.train.prelim)
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

glob <- gtKEGG( Surv(exp(pheno.train[,3]), pheno.train[,2]), exprs.norm, annotation = 'org.Hs.eg.db', multtest = "BH",probe2entrez = xx)


pathway.names <- names(glob)[1:30]


yy <- as.list(KEGGPATHID2EXTID)

######### TO GET GENE NAMES of GENES WITHIN PATHWAYS ##############################3

# For the reverse map:
# Convert the object to a list
dd <- as.list(org.Hs.egPATH2EG)
# Remove pathway identifiers that do not map to any entrez gene id
dd <- dd[!is.na(dd)]
if(length(dd) > 0){
  # The entrez gene identifiers for the first two elements of XX
  dd[1:2]
  # Get the first one
  dd[[1]]
}

path2entrez <- dd

path2entrez.subset <- path2entrez[pathway.names]

library('biomaRt')
# listDatasets(ensembl)
ensembl = useMart('ensembl')
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

names.genes <- list(0)
for( i in 1:length(pathway.names)){
  names.genes[[i]] = getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "entrezgene",values =path2entrez.subset[i], mart=ensembl)
}

### Data frame FROM the EXPRESSION #################

data.exprs <- as.data.frame(exprs.norm)
data.subsets <- list(0)
for ( i in 1:length(pathway.names) ){
  
  temp.names <- names.genes[[i]][,2][c(names.genes[[i]][,2])%in% colnames(data.exprs)]
  data.subsets[[i]] <- as.data.frame(data.exprs[, temp.names]) 
}

data.combined <- matrix(0, nrow = nrow(data.exprs), ncol =length(pathway.names) )
for ( i in 1:length(pathway.names)){
  #   data.combined[,i] <- as.vector(apply(data.subsets[[i]],1, mean))
  pc <- prcomp(data.subsets[[i]])
  data.combined[,i] <- predict(pc,newdata = data.subsets[[i]])[,1]
}

rownames(data.combined) <- rownames(data.exprs)
colnames(data.combined) <- pathway.names

############ A Small Vizualization ###################
pc <- prcomp(data.combined)
pc.pred <- predict(pc,newdata = data.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = pheno.train[,4])


### checking the Heatmap #######################
### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = t(data.combined), scale = "row",Rowv =TRUE ,Colv = FALSE, dendrogram = "row", ColSideColors =c("blue","red","green","purple")[c.true], labRow = colnames(data.combined), labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")



######### Saving the Top 30 pathways #################
save(data.combined,file = '30pathwayVerhaark.RData')

