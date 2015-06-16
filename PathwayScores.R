############## This Code Calculates the Pathway Scores for the different data sets #############
############ VERHAARK DATA SET #############################################################
setwd('/home/bit/ashar/ExpressionSets/Verhark')
exprs <- read.csv('/home/bit/ashar/ExpressionSets/Verhark/UNC202.txt', sep = '\t', header = FALSE)
patient.names <- exprs[1,]
exprs <- exprs[-1,]
patient.names <- patient.names[-1]
names <- as.character(unlist(patient.names))

#### Removing some dirty names and replacing with good names###########
q <- strsplit(names,"-01A")
names.q <- c(0)
for ( i in 1:length(q)){
  names.q <- c(names.q,q[[i]][1]) 
}
names.q <- gsub("-01B-01","",names.q)
names.q <- gsub("-01C-01","",names.q)
names.q <- names.q[-1]


## Read PhenoData
library('xlsx')
pheno <- read.xlsx("survival.xls", sheetIndex =1)
pheno[,1] <- as.character(as.matrix(pheno[,1]))


## Arranging the expression object to have the same order as the 
exprs.ready <- exprs[2:11862,2:203]
rownames(exprs.ready) <- as.character(exprs[2:11862,1])
colnames(exprs.ready) <- names.q

exprs.norm <- t(exprs.ready)

mt <- match(rownames(exprs.norm),pheno[,1])
pheno <- pheno[mt,] 


## Taking the log transform of the data
ind.na <- which(as.vector(pheno[,3])== "NA")
ind.nna <- which(as.vector(pheno[,3])!= "NA")

templog <- as.numeric(as.vector(pheno[ind.nna,3]))
templog <- log(templog)

pheno[,3] <- as.vector(pheno[,3])

pheno[ind.nna,3] <- templog
pheno[ind.na,3] <- NA

pheno[,3] <- as.numeric(as.matrix(pheno[,3]))

## Expression Object is Ready


## Delete Those Patients which have NA
pheno <- pheno[-ind.na,]
exprs.norm <- exprs.norm[-ind.na,]


surv <- Surv(time = pheno[,3], event = pheno[,2])


library(globaltest)
library(org.Hs.eg.db)
library(GO.db)
library(KEGG.db)
exprs.norm <- as.matrix(as.data.frame(exprs.norm))
class(exprs.norm) <- "numeric"
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

glob <- gtKEGG( pheno[,3], (exprs.norm), annotation = 'org.Hs.eg.db', multtest = "BH",probe2entrez = xx)


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
pdf("VizVerhaarkDataPCA.pdf")
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = pheno[,4])
dev.off()


########################################################################
#########################################################################

save(data.combined,file = '30pathway.RData')
