---
title: "Minute Zika virus comparison"
date: "May 29, 2019"
output:
  html_document: 
    fig_height: 25
    fig_width: 25
    highlight: tango
    keep_md: yes
    theme: cosmo
    toc: yes
---

##Minute Zika virus comparison
Analysis of gene expression pattern of ZikaV data. Normal patients vs Zika Virus using weighted gene co-expression network analysis (WGCNA).



  
Packages used in this study:

```r
library(tidyverse)
library(knitr)
library(WGCNA)
library(flashClust)
library(WGCNA)
library(psycho)
library(tximport)
library(tibble)
library(ggplot2)
library(tximport)
library(kableExtra)
library(ggraph)
library(igraph)
library(tidyverse)
library(IRanges)
library(clusterProfiler)
library(viridis)
library(S4Vectors)
library(circlize)
library(chorddiag)
options(stringsAsFactors = FALSE)# Important for WGCNA!
```


#### 1- Data input and cleaning
Load metadata information

##### Minute Zika virus comparison


```r
Metadata <- read_csv("Human_genome/METADATA.csv")
```

```
## Parsed with column specification:
## cols(
##   Samples = col_character(),
##   Treatment = col_character(),
##   tissue = col_character(),
##   Accession = col_character(),
##   BioProject = col_character(),
##   Study = col_character()
## )
```

```r
#head(Metadata)

#Parcel metadata file to use only data from BioProject PRJNA497590
PRJNA497590_metadata <- Metadata %>% filter(BioProject %in% c("PRJNA497590")) %>% group_by(Samples) 

#Set working directories for the reads
PRJNA497590_files <- file.path("Human_genome","PRJNA497590", paste0(PRJNA497590_metadata$Samples, ".genes.results"))
all(file.exists(PRJNA497590_files))
```

```
## [1] TRUE
```

```r
#Rename reads (matching information on the metadata)
names(PRJNA497590_files ) <- c(PRJNA497590_metadata$Samples)

PRJNA497590_count <- tximport(PRJNA497590_files, type = "rsem", txOut = TRUE, txIdCol = "Name", abundanceCol = "TPM") # to use in WGCNA
```

```
## It looks like you are importing RSEM genes.results files, setting txIn=FALSE
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5
```

```r
#Isolate RSEM counts

PRJNA497590_counts<-as.data.frame(PRJNA497590_count$counts)%>%
  rownames_to_column() %>% 
   dplyr::select(human_gene_id = rowname, everything())
#head(PRJNA497590_counts)

#Save counts
write.csv(PRJNA497590_counts, "Human_genome/PRJNA497590_counts.csv")
```
##### Minute data set

```r
Fly_treatments <- read_csv("Fly/METADATA/detailed_samples.csv", col_types = cols()) %>% replace(.=="a_w", "wild_type") %>% replace(.=="w_", "control")

#Set working directories for the reads
Files <- file.path("Fly", "RSEM", paste0(Fly_treatments$sample, ".genes.results"))

#Rename reads (matching information on the metadata)
names(Files) <- c(Fly_treatments$sample)
Gene_counts <- tximport(Files, type = "rsem", txOut = TRUE, txIdCol = "Name", abundanceCol = "TPM") # to use in WGCNA

#Isolate RSEM counts

Gene_count <- as.data.frame(Gene_counts$counts)
#head(Gene_count)
```
Filter all genes that were differentially expressed between wild type and control 

```r
dge_list <-read_csv("Wild_sig.csv") %>% 
  rename("target_id" = "gene_id") %>%
  filter(padj < 0.05) #making sure all p-adjusted value are smaller than 0.05
```

```
## Parsed with column specification:
## cols(
##   target_id = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double(),
##   ext_gene = col_character()
## )
```

```r
'%ni%' <- Negate('%in%') ##function for excluding wildtype results

Gene_count <- Gene_count %>% 
  filter(rownames(Gene_count) %ni% dge_list$gene_id)

#Save counts
write.csv(Gene_count, "Fly/Gene_count.csv")
```

#### 2- Isolate all orthologus genes between _Homo sapiens_ and _Drosophila melanogaster_ 
Load ensembl data set

```r
# ftp://ftp.flybase.org/releases/FB2019_02/precomputed_files/orthologs/dmel_human_orthologs_disease_fb_2019_02.tsv.gz
human_ortholog <- read_csv("Human_genome/dmel_human_orthologs.csv", col_types = cols()) ##Downloaded from fly base 
human_id <- read_csv("Human_genome/human_id.csv", col_types = cols()) ##Downloaded from ensembl
Dml_Human_orthologs <- human_id %>% dplyr::left_join(human_ortholog,by= "HGNC_ID") %>% dplyr::select(fly_gene_ID, human_gene_id = gene_id, fly_gene_symbol) %>% drop
```
##### BioProject: PRJNA497590

```r
#Isolate orthologus genes between _Homo sapiens_ and _Drosophila melanogaster_
PRJNA497590_counts <- PRJNA497590_counts %>% dplyr::left_join(Dml_Human_orthologs, by = "human_gene_id") %>% drop_na(fly_gene_ID) %>% distinct(fly_gene_ID, .keep_all = TRUE) %>% dplyr::select(fly_gene_ID, c(PRJNA497590_metadata$Samples))
```
##### Minute project

```r
Minute_counts <- read_csv("Fly/Gene_count.csv") %>% rename_(fly_gene_ID = names(.)[1])
#Isolate orthologus genes between _Homo sapiens_ and _Drosophila melanogaster_
#Minute_counts <- Gene_count %>% left_join(Dml_Human_orthologs, by = "fly_gene_ID") %>% drop_na(fly_gene_ID) %>% distinct(fly_gene_ID, .keep_all = TRUE) %>% dplyr::select(fly_gene_ID, c(PRJNA497590_metadata$Samples))
```
Loading the gene names and descriptions

```r
Fly_g0 <- read_tsv("Fly/METADATA/mart_export.txt") %>% dplyr::select(Fly_go = GO_term_name, gene_id,gene_name) %>% drop_na()
```

#### 3- Pre-processing data sets (ensures compatibility)
Isolate common genes

```r
common_gene_id <- unique(Reduce(intersect,list(PRJNA497590_counts$fly_gene_ID, Minute_counts$fly_gene_ID)))
PRJNA497590_counts = PRJNA497590_counts %>% filter(fly_gene_ID %in% common_gene_id) %>% group_by(fly_gene_ID) %>% arrange(fly_gene_ID) 
Minute_counts = Minute_counts %>% filter(fly_gene_ID %in% common_gene_id) %>% dplyr::select(fly_gene_ID, c(Fly_treatments$sample)) %>% group_by(fly_gene_ID) %>% arrange(fly_gene_ID)
```

## Minute and Zika Virus comparison  {.tabset .tabset-fade}

### The _RpL3_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPL3_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPL3", "control"))
RPL3_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPL3_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpL3")

RPL3_gene_expression = vector(mode = "list", length = nSets)

RPL3_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPL3_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPL3_gene_expression[[2]] = list(data = as.data.frame(t(RPL3_counts[,-1])));
names(RPL3_gene_expression[[2]]$data) = RPL3_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPL3_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPL3_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
  {# Remove the offending genes and samples
    RPL3_gene_expression[[set]]$data = RPL3_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  } # Update exprSize
  exprSize = checkSets(RPL3_gene_expression)
  gsg2 = goodSamplesGenesMS(RPL3_gene_expression, verbose = 2); #check again
  gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```
Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPL3_gene_expression)
RPL3_softpower = softpower

RPL3_human_net<-network(RPL3_gene_expression[[1]]$data)
RPL3_fly_net<- network(RPL3_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)

```r
# Human network:
RPL3_AdjMatHuman =  abs(cor(RPL3_gene_expression[[1]]$data,use="p"))^RPL3_softpower

# Fly network:
RPL3_AdjMatFly =  abs(cor(RPL3_gene_expression[[2]]$data,use="p"))^RPL3_softpower

## Calculation of the whole network connectivity k:
RPL3_ConnectivityHuman=as.vector(apply(RPL3_AdjMatHuman,2,sum, na.rm=T))

RPL3_ConnectivityFly =as.vector(apply(RPL3_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPL3_ConnectivityHuman=RPL3_ConnectivityHuman %>% psycho::standardize() 
RPL3_ConnectivityFly=RPL3_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPL3_ExpressionHuman = apply(RPL3_gene_expression[[1]]$data,2,mean)
RPL3_ExpressionHuman = RPL3_ExpressionHuman/max(RPL3_ExpressionHuman)
RPL3_ExpressionFly = apply(RPL3_gene_expression[[2]]$data,2,mean)
RPL3_ExpressionFly= RPL3_ExpressionFly/max(RPL3_ExpressionFly)
RPL3_cor_data<- data.frame(cbind(RPL3_ExpressionHuman, RPL3_ExpressionFly))
                      
cor.test(RPL3_cor_data$RPL3_ExpressionHuman,RPL3_cor_data$RPL3_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL3_cor_data$RPL3_ExpressionHuman and RPL3_cor_data$RPL3_ExpressionFly
## S = 33909063291, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3078535
```

```r
ggplot(RPL3_cor_data, aes(RPL3_ExpressionHuman, RPL3_ExpressionFly)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPL3_ConnectivityHuman,RPL3_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL3_ConnectivityHuman and RPL3_ConnectivityFly
## S = 48648534425, p-value = 0.5686
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##         rho 
## 0.006993657
```

```r
RPL3_net_data<- data.frame(cbind(RPL3_ConnectivityHuman,RPL3_ConnectivityFly))

ggplot(RPL3_net_data, aes(RPL3_ConnectivityHuman, RPL3_ConnectivityFly)) + 
  geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
  geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpL3 in the Minute data

```r
#Isolate RpL3 gene module
RPL3_Module <- RPL3_human_net$colors[which(colnames(RPL3_gene_expression[[1]]$data) == "FBgn0020910")]
RPL3_ModuleGenes <- colnames(RPL3_gene_expression[[1]]$data)[RPL3_human_net$colors == RPL3_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaVirus", "Experiment") 
ZikaV_RPL3_expression  <- list(Experiment = list(data = RPL3_gene_expression[[2]]$data), ZikaV = list(data = RPL3_gene_expression[[1]]$data))
ZikaV_RPL3_multiColor <- list(Experiment = labels2colors(RPL3_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPL3_mp<- modulePreservation(ZikaV_RPL3_expression , ZikaV_RPL3_multiColor, nPermutations = npermutation,
                                   referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPL3_mp, file = "RDS/ZikaV_RPL3_mp.RData");

#_P-value_

ZikaV_RPL3_p.value=exp(ZikaV_RPL3_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL3_Module),"log.psummary.pres"])
ZikaV_RPL3_Zsummary = ZikaV_RPL3_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL3_Module),"Zsummary.pres"]

ZikaV_mp= cbind(Zsummary=ZikaV_RPL3_Zsummary, p.value=ZikaV_RPL3_p.value)
rownames(ZikaV_mp)[1] <- "RpL3"
```
Comparing correlation of the RpL3 module between humans and the Minute project


```r
cor.test(RPL3_ExpressionHuman[RPL3_ModuleGenes],RPL3_ExpressionFly[RPL3_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL3_ExpressionHuman[RPL3_ModuleGenes] and RPL3_ExpressionFly[RPL3_ModuleGenes]
## S = 923450, p-value = 0.000000001513
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.4017032
```

```r
RPL3_shared_module<- data.frame(cbind(RPL3_ExpressionHuman = RPL3_ExpressionHuman[RPL3_ModuleGenes], RPL3_ExpressionFly= RPL3_ExpressionFly[RPL3_ModuleGenes]))

ggplot(RPL3_shared_module, aes(RPL3_ExpressionHuman, RPL3_ExpressionFly)) + 
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

### The _RpL14_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPL14_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPL14", "control"))
RPL14_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPL14_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpL14")

RPL14_gene_expression = vector(mode = "list", length = nSets)

RPL14_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPL14_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPL14_gene_expression[[2]] = list(data = as.data.frame(t(RPL14_counts[,-1])));
names(RPL14_gene_expression[[2]]$data) = RPL14_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPL14_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPL14_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPL14_gene_expression[[set]]$data = RPL14_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPL14_gene_expression)
   gsg2 = goodSamplesGenesMS(RPL14_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```
Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPL14_gene_expression)
RPL14_softpower = softpower

RPL14_human_net<-network(RPL14_gene_expression[[1]]$data)
RPL14_fly_net<- network(RPL14_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)


```r
# Human network:
RPL14_AdjMatHuman =  abs(cor(RPL14_gene_expression[[1]]$data,use="p"))^RPL14_softpower

# Fly network:
RPL14_AdjMatFly =  abs(cor(RPL14_gene_expression[[2]]$data,use="p"))^RPL14_softpower

## Calculation of the whole network connectivity k:
RPL14_ConnectivityHuman=as.vector(apply(RPL14_AdjMatHuman,2,sum, na.rm=T))

RPL14_ConnectivityFly =as.vector(apply(RPL14_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPL14_ConnectivityHuman=RPL14_ConnectivityHuman %>% psycho::standardize() 
RPL14_ConnectivityFly=RPL14_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPL14_ExpressionHuman = apply(RPL14_gene_expression[[1]]$data,2,mean)
RPL14_ExpressionHuman = RPL14_ExpressionHuman/max(RPL14_ExpressionHuman)
RPL14_ExpressionFly = apply(RPL14_gene_expression[[2]]$data,2,mean)
RPL14_ExpressionFly= RPL14_ExpressionFly/max(RPL14_ExpressionFly)
RPL14_cor_data<- data.frame(cbind(RPL14_ExpressionHuman, RPL14_ExpressionFly))

cor.test(RPL14_cor_data$RPL14_ExpressionHuman,RPL14_cor_data$RPL14_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL14_cor_data$RPL14_ExpressionHuman and RPL14_cor_data$RPL14_ExpressionFly
## S = 34072331429, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##      rho 
## 0.301689
```

```r
ggplot(RPL14_cor_data, aes(RPL14_ExpressionHuman, RPL14_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPL14_ConnectivityHuman,RPL14_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL14_ConnectivityHuman and RPL14_ConnectivityFly
## S = 49113024805, p-value = 0.5925
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##          rho 
## -0.006569356
```

```r
RPL14_net_data<- data.frame(cbind(RPL14_ConnectivityHuman,RPL14_ConnectivityFly))

ggplot(RPL14_net_data, aes(RPL14_ConnectivityHuman, RPL14_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpL14 in the Minute data

```r
#Isolate RpL14 gene module
RPL14_Module <- RPL14_human_net$colors[which(colnames(RPL14_gene_expression[[1]]$data) == "FBgn0017579")]
RPL14_ModuleGenes <- colnames(RPL14_gene_expression[[1]]$data)[RPL14_human_net$colors == RPL14_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPL14_expression  <- list(Experiment = list(data = RPL14_gene_expression[[2]]$data), ZikaV = list(data = RPL14_gene_expression[[1]]$data))
ZikaV_RPL14_multiColor <- list(Experiment = labels2colors(RPL14_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPL14_mp<- modulePreservation(ZikaV_RPL14_expression , ZikaV_RPL14_multiColor, nPermutations = npermutation,
                                          referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPL14_mp, file = "RDS/ZikaV_RPL14_mp.RData");

#_P-value_


ZikaV_RPL14_p.value=exp(ZikaV_RPL14_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL14_Module),"log.psummary.pres"])
ZikaV_RPL14_Zsummary = ZikaV_RPL14_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL14_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpL14" = c(ZikaV_RPL14_Zsummary,ZikaV_RPL14_p.value))
```
Comparing correlation of the RpL14 module between humans and the Minute project


```r
cor.test(RPL14_ExpressionHuman[RPL14_ModuleGenes],RPL14_ExpressionFly[RPL14_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL14_ExpressionHuman[RPL14_ModuleGenes] and RPL14_ExpressionFly[RPL14_ModuleGenes]
## S = 70911079, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3216643
```

```r
RPL14_shared_module<- data.frame(cbind(RPL14_ExpressionHuman = RPL14_ExpressionHuman[RPL14_ModuleGenes], RPL14_ExpressionFly= RPL14_ExpressionFly[RPL14_ModuleGenes]))

ggplot(RPL14_shared_module, aes(RPL14_ExpressionHuman, RPL14_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

### The _RpS13_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPS13_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPS13", "control"))
RPS13_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPS13_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpS13")

RPS13_gene_expression = vector(mode = "list", length = nSets)

RPS13_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPS13_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPS13_gene_expression[[2]] = list(data = as.data.frame(t(RPS13_counts[,-1])));
names(RPS13_gene_expression[[2]]$data) = RPS13_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPS13_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPS13_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPS13_gene_expression[[set]]$data = RPS13_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPS13_gene_expression)
   gsg2 = goodSamplesGenesMS(RPS13_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```

Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPS13_gene_expression)
RPS13_softpower = softpower

RPS13_human_net<-network(RPS13_gene_expression[[1]]$data)
RPS13_fly_net<- network(RPS13_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)

```r
# Human network:
 RPS13_AdjMatHuman =  abs(cor( RPS13_gene_expression[[1]]$data,use="p"))^ RPS13_softpower

# Fly network:
 RPS13_AdjMatFly =  abs(cor( RPS13_gene_expression[[2]]$data,use="p"))^ RPS13_softpower

## Calculation of the whole network connectivity k:
 RPS13_ConnectivityHuman=as.vector(apply( RPS13_AdjMatHuman,2,sum, na.rm=T))

 RPS13_ConnectivityFly =as.vector(apply( RPS13_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
 RPS13_ConnectivityHuman= RPS13_ConnectivityHuman %>% psycho::standardize() 
 RPS13_ConnectivityFly= RPS13_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
 RPS13_ExpressionHuman = apply( RPS13_gene_expression[[1]]$data,2,mean)
 RPS13_ExpressionHuman =  RPS13_ExpressionHuman/max( RPS13_ExpressionHuman)
 RPS13_ExpressionFly = apply(RPS13_gene_expression[[2]]$data,2,mean)
 RPS13_ExpressionFly=  RPS13_ExpressionFly/max( RPS13_ExpressionFly)
 RPS13_cor_data<- data.frame(cbind( RPS13_ExpressionHuman,  RPS13_ExpressionFly))

cor.test( RPS13_cor_data$ RPS13_ExpressionHuman, RPS13_cor_data$ RPS13_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS13_cor_data$RPS13_ExpressionHuman and RPS13_cor_data$RPS13_ExpressionFly
## S = 34442590088, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##      rho 
## 0.299179
```

```r
ggplot( RPS13_cor_data, aes( RPS13_ExpressionHuman,  RPS13_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test( RPS13_ConnectivityHuman, RPS13_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS13_ConnectivityHuman and RPS13_ConnectivityFly
## S = 47276403101, p-value = 0.001908
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.0380428
```

```r
 RPS13_net_data<- data.frame(cbind( RPS13_ConnectivityHuman, RPS13_ConnectivityFly))

ggplot( RPS13_net_data, aes( RPS13_ConnectivityHuman,  RPS13_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-40-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpS13 in the Minute data

```r
#Isolate RpS13 gene module
RPS13_Module <- RPS13_human_net$colors[which(colnames(RPS13_gene_expression[[1]]$data) == "FBgn0010265")]
RPS13_ModuleGenes <- colnames(RPS13_gene_expression[[1]]$data)[RPS13_human_net$colors == RPS13_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPS13_expression  <- list(Experiment = list(data = RPS13_gene_expression[[2]]$data), ZikaV = list(data = RPS13_gene_expression[[1]]$data))
ZikaV_RPS13_multiColor <- list(Experiment = labels2colors(RPS13_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPS13_mp<- modulePreservation(ZikaV_RPS13_expression , ZikaV_RPS13_multiColor, nPermutations = npermutation,
                                          referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPS13_mp, file = "RDS/ZikaV_RPS13_mp.RData");

#_P-value_
ZikaV_RPS13_p.value=exp(ZikaV_RPS13_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPS13_Module),"log.psummary.pres"])
ZikaV_RPS13_Zsummary = ZikaV_RPS13_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPS13_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpS13" = c(ZikaV_RPS13_Zsummary,ZikaV_RPS13_p.value))
```
Comparing correlation of the RpS13 module between humans and the Minute project


```r
cor.test(RPS13_ExpressionHuman[RPS13_ModuleGenes],RPS13_ExpressionFly[RPS13_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS13_ExpressionHuman[RPS13_ModuleGenes] and RPS13_ExpressionFly[RPS13_ModuleGenes]
## S = 2928756606, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.2400986
```

```r
RPS13_shared_module<- data.frame(cbind(RPS13_ExpressionHuman = RPS13_ExpressionHuman[RPS13_ModuleGenes], RPS13_ExpressionFly= RPS13_ExpressionFly[RPS13_ModuleGenes]))

ggplot(RPS13_shared_module, aes(RPS13_ExpressionHuman, RPS13_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-44-1.png" style="display: block; margin: auto;" />

### The _RpL19_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPL19_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPL19", "control"))
RPL19_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPL19_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpL19")

RPL19_gene_expression = vector(mode = "list", length = nSets)

RPL19_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPL19_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPL19_gene_expression[[2]] = list(data = as.data.frame(t(RPL19_counts[,-1])));
names(RPL19_gene_expression[[2]]$data) = RPL19_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPL19_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPL19_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPL19_gene_expression[[set]]$data = RPL19_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPL19_gene_expression)
   gsg2 = goodSamplesGenesMS(RPL19_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```

Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPL19_gene_expression)
RPL19_softpower = softpower

RPL19_human_net<-network(RPL19_gene_expression[[1]]$data)
RPL19_fly_net<- network(RPL19_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)

```r
# Human network:
RPL19_AdjMatHuman =  abs(cor(RPL19_gene_expression[[1]]$data,use="p"))^RPL19_softpower

# Fly network:
RPL19_AdjMatFly =  abs(cor(RPL19_gene_expression[[2]]$data,use="p"))^RPL19_softpower

## Calculation of the whole network connectivity k:
RPL19_ConnectivityHuman=as.vector(apply(RPL19_AdjMatHuman,2,sum, na.rm=T))

RPL19_ConnectivityFly =as.vector(apply(RPL19_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPL19_ConnectivityHuman=RPL19_ConnectivityHuman %>% psycho::standardize() 
RPL19_ConnectivityFly=RPL19_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPL19_ExpressionHuman = apply(RPL19_gene_expression[[1]]$data,2,mean)
RPL19_ExpressionHuman = RPL19_ExpressionHuman/max(RPL19_ExpressionHuman)
RPL19_ExpressionFly = apply(RPL19_gene_expression[[2]]$data,2,mean)
RPL19_ExpressionFly= RPL19_ExpressionFly/max(RPL19_ExpressionFly)
RPL19_cor_data<- data.frame(cbind(RPL19_ExpressionHuman, RPL19_ExpressionFly))

cor.test(RPL19_cor_data$RPL19_ExpressionHuman,RPL19_cor_data$RPL19_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL19_cor_data$RPL19_ExpressionHuman and RPL19_cor_data$RPL19_ExpressionFly
## S = 34655163296, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.2983382
```

```r
ggplot(RPL19_cor_data, aes(RPL19_ExpressionHuman, RPL19_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-51-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPL19_ConnectivityHuman,RPL19_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL19_ConnectivityHuman and RPL19_ConnectivityFly
## S = 49437931127, p-value = 0.937
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##           rho 
## -0.0009679751
```

```r
RPL19_net_data<- data.frame(cbind(RPL19_ConnectivityHuman,RPL19_ConnectivityFly))

ggplot(RPL19_net_data, aes(RPL19_ConnectivityHuman, RPL19_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-52-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpL19 in the Minute data

```r
#Isolate RpL19 gene module
RPL19_Module <- RPL19_human_net$colors[which(colnames(RPL19_gene_expression[[1]]$data) == "FBgn0002607")]
RPL19_ModuleGenes <- colnames(RPL19_gene_expression[[1]]$data)[RPL19_human_net$colors == RPL19_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPL19_expression  <- list(Experiment = list(data = RPL19_gene_expression[[2]]$data), ZikaV = list(data = RPL19_gene_expression[[1]]$data))
ZikaV_RPL19_multiColor <- list(Experiment = labels2colors(RPL19_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPL19_mp<- modulePreservation(ZikaV_RPL19_expression , ZikaV_RPL19_multiColor, nPermutations = npermutation,
                                            referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPL19_mp, file = "RDS/ZikaV_RPL19_mp.RData");

#_P-value_
ZikaV_RPL19_p.value=exp(ZikaV_RPL19_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL19_Module),"log.psummary.pres"])
ZikaV_RPL19_Zsummary = ZikaV_RPL19_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL19_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpL19" = c(ZikaV_RPL19_Zsummary,ZikaV_RPL19_p.value))
```
Comparing correlation of the RpL19 module between humans and the Minute project


```r
cor.test(RPL19_ExpressionHuman[RPL19_ModuleGenes],RPL19_ExpressionFly[RPL19_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL19_ExpressionHuman[RPL19_ModuleGenes] and RPL19_ExpressionFly[RPL19_ModuleGenes]
## S = 345465, p-value = 0.000006739
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3605744
```

```r
RPL19_shared_module<- data.frame(cbind(RPL19_ExpressionHuman = RPL19_ExpressionHuman[RPL19_ModuleGenes], RPL19_ExpressionFly= RPL19_ExpressionFly[RPL19_ModuleGenes]))

ggplot(RPL19_shared_module, aes(RPL19_ExpressionHuman, RPL19_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-56-1.png" style="display: block; margin: auto;" />

### The _RpS19b_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPS19b_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPS19b", "control"))
RPS19b_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPS19b_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpS19b")

RPS19b_gene_expression = vector(mode = "list", length = nSets)

RPS19b_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPS19b_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPS19b_gene_expression[[2]] = list(data = as.data.frame(t(RPS19b_counts[,-1])));
names(RPS19b_gene_expression[[2]]$data) = RPS19b_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPS19b_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPS19b_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPS19b_gene_expression[[set]]$data = RPS19b_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPS19b_gene_expression)
   gsg2 = goodSamplesGenesMS(RPS19b_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```

Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPS19b_gene_expression)
RPS19b_softpower = softpower

RPS19b_human_net<-network(RPS19b_gene_expression[[1]]$data)
RPS19b_fly_net<- network(RPS19b_gene_expression[[2]]$data)
```

#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)


```r
# Human network:
RPS19b_AdjMatHuman =  abs(cor(RPS19b_gene_expression[[1]]$data,use="p"))^RPS19b_softpower

# Fly network:
RPS19b_AdjMatFly =  abs(cor(RPS19b_gene_expression[[2]]$data,use="p"))^RPS19b_softpower

## Calculation of the whole network connectivity k:
RPS19b_ConnectivityHuman=as.vector(apply(RPS19b_AdjMatHuman,2,sum, na.rm=T))

RPS19b_ConnectivityFly =as.vector(apply(RPS19b_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPS19b_ConnectivityHuman=RPS19b_ConnectivityHuman %>% psycho::standardize() 
RPS19b_ConnectivityFly=RPS19b_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPS19b_ExpressionHuman = apply(RPS19b_gene_expression[[1]]$data,2,mean)
RPS19b_ExpressionHuman = RPS19b_ExpressionHuman/max(RPS19b_ExpressionHuman)
RPS19b_ExpressionFly = apply(RPS19b_gene_expression[[2]]$data,2,mean)
RPS19b_ExpressionFly= RPS19b_ExpressionFly/max(RPS19b_ExpressionFly)
RPS19b_cor_data<- data.frame(cbind(RPS19b_ExpressionHuman, RPS19b_ExpressionFly))

cor.test(RPS19b_cor_data$RPS19b_ExpressionHuman,RPS19b_cor_data$RPS19b_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS19b_cor_data$RPS19b_ExpressionHuman and RPS19b_cor_data$RPS19b_ExpressionFly
## S = 33837427706, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3077552
```

```r
ggplot(RPS19b_cor_data, aes(RPS19b_ExpressionHuman, RPS19b_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-63-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPS19b_ConnectivityHuman,RPS19b_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS19b_ConnectivityHuman and RPS19b_ConnectivityFly
## S = 48840886123, p-value = 0.947
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##          rho 
## 0.0008149587
```

```r
RPS19b_net_data<- data.frame(cbind(RPS19b_ConnectivityHuman,RPS19b_ConnectivityFly))

ggplot(RPS19b_net_data, aes(RPS19b_ConnectivityHuman, RPS19b_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-64-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpS19b in the Minute data

```r
#Isolate RpS19b gene module
RPS19b_Module <- RPS19b_human_net$colors[which(colnames(RPS19b_gene_expression[[1]]$data) == "FBgn0039129")]
RPS19b_ModuleGenes <- colnames(RPS19b_gene_expression[[1]]$data)[RPS19b_human_net$colors == RPS19b_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPS19b_expression  <- list(Experiment = list(data = RPS19b_gene_expression[[2]]$data), ZikaV = list(data = RPS19b_gene_expression[[1]]$data))
ZikaV_RPS19b_multiColor <- list(Experiment = labels2colors(RPS19b_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPS19b_mp<- modulePreservation(ZikaV_RPS19b_expression , ZikaV_RPS19b_multiColor, nPermutations = npermutation,
                                          referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPS19b_mp, file = "RDS/ZikaV_RPS19b_mp.RData");

#_P-value_
ZikaV_RPS19b_p.value=exp(ZikaV_RPS19b_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPS19b_Module),"log.psummary.pres"])
ZikaV_RPS19b_Zsummary = ZikaV_RPS19b_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPS19b_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpS19b" = c(ZikaV_RPS19b_Zsummary,ZikaV_RPS19b_p.value))
```
Comparing correlation of the RpS19b module between humans and the Minute project


```r
cor.test(RPS19b_ExpressionHuman[RPS19b_ModuleGenes],RPS19b_ExpressionFly[RPS19b_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPS19b_ExpressionHuman[RPS19b_ModuleGenes] and RPS19b_ExpressionFly[RPS19b_ModuleGenes]
## S = 7990533447, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.2707548
```

```r
RPS19b_shared_module<- data.frame(cbind(RPS19b_ExpressionHuman = RPS19b_ExpressionHuman[RPS19b_ModuleGenes], RPS19b_ExpressionFly= RPS19b_ExpressionFly[RPS19b_ModuleGenes]))

ggplot(RPS19b_shared_module, aes(RPS19b_ExpressionHuman, RPS19b_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-68-1.png" style="display: block; margin: auto;" />
### The _RpL24-like _- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPL24_like_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPL24_like ", "control"))
RPL24_like_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPL24_like_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpL24-like ")

RPL24_like_gene_expression = vector(mode = "list", length = nSets)

RPL24_like_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPL24_like_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPL24_like_gene_expression[[2]] = list(data = as.data.frame(t(RPL24_like_counts[,-1])));
names(RPL24_like_gene_expression[[2]]$data) = RPL24_like_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPL24_like_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPL24_like_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPL24_like_gene_expression[[set]]$data = RPL24_like_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPL24_like_gene_expression)
   gsg2 = goodSamplesGenesMS(RPL24_like_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```

Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPL24_like_gene_expression)
RPL24_like_softpower = softpower

RPL24_like_human_net<-network(RPL24_like_gene_expression[[1]]$data)
RPL24_like_fly_net<- network(RPL24_like_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)

```r
# Human network:
RPL24_like_AdjMatHuman =  abs(cor(RPL24_like_gene_expression[[1]]$data,use="p"))^RPL24_like_softpower

# Fly network:
RPL24_like_AdjMatFly =  abs(cor(RPL24_like_gene_expression[[2]]$data,use="p"))^RPL24_like_softpower

## Calculation of the whole network connectivity k:
RPL24_like_ConnectivityHuman=as.vector(apply(RPL24_like_AdjMatHuman,2,sum, na.rm=T))

RPL24_like_ConnectivityFly =as.vector(apply(RPL24_like_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPL24_like_ConnectivityHuman=RPL24_like_ConnectivityHuman %>% psycho::standardize() 
RPL24_like_ConnectivityFly=RPL24_like_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPL24_like_ExpressionHuman = apply(RPL24_like_gene_expression[[1]]$data,2,mean)
RPL24_like_ExpressionHuman = RPL24_like_ExpressionHuman/max(RPL24_like_ExpressionHuman)
RPL24_like_ExpressionFly = apply(RPL24_like_gene_expression[[2]]$data,2,mean)
RPL24_like_ExpressionFly= RPL24_like_ExpressionFly/max(RPL24_like_ExpressionFly)
RPL24_like_cor_data<- data.frame(cbind(RPL24_like_ExpressionHuman, RPL24_like_ExpressionFly))

cor.test(RPL24_like_cor_data$RPL24_like_ExpressionHuman,RPL24_like_cor_data$RPL24_like_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL24_like_cor_data$RPL24_like_ExpressionHuman and RPL24_like_cor_data$RPL24_like_ExpressionFly
## S = 32582887813, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3052632
```

```r
ggplot(RPL24_like_cor_data, aes(RPL24_like_ExpressionHuman, RPL24_like_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-75-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPL24_like_ConnectivityHuman,RPL24_like_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL24_like_ConnectivityHuman and RPL24_like_ConnectivityFly
## S = 47466892844, p-value = 0.3276
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##         rho 
## -0.01209566
```

```r
RPL24_like_net_data<- data.frame(cbind(RPL24_like_ConnectivityHuman,RPL24_like_ConnectivityFly))

ggplot(RPL24_like_net_data, aes(RPL24_like_ConnectivityHuman, RPL24_like_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-76-1.png" style="display: block; margin: auto;" />
Identify the genes of interest on module RpL24-like  in the Minute data

```r
#Isolate RpL24-like  gene module
RPL24_like_Module <- RPL24_like_human_net$colors[which(colnames(RPL24_like_gene_expression[[1]]$data) == "FBgn0037899")]
RPL24_like_ModuleGenes <- colnames(RPL24_like_gene_expression[[1]]$data)[RPL24_like_human_net$colors == RPL24_like_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPL24_like_expression  <- list(Experiment = list(data = RPL24_like_gene_expression[[2]]$data), ZikaV = list(data = RPL24_like_gene_expression[[1]]$data))
ZikaV_RPL24_like_multiColor <- list(Experiment = labels2colors(RPL24_like_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPL24_like_mp<- modulePreservation(ZikaV_RPL24_like_expression , ZikaV_RPL24_like_multiColor, nPermutations = npermutation,
                                          referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPL24_like_mp, file = "RDS/ZikaV_RPL24_like_mp.RData");

#_P-value_

ZikaV_RPL24_like_p.value=exp(ZikaV_RPL24_like_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL24_like_Module),"log.psummary.pres"])
ZikaV_RPL24_like_Zsummary = ZikaV_RPL24_like_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL24_like_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpL24-like" = c(ZikaV_RPL24_like_Zsummary,ZikaV_RPL24_like_p.value))
```
Comparing correlation of the RpL24-like  module between humans and the Minute project


```r
cor.test(RPL24_like_ExpressionHuman[RPL24_like_ModuleGenes],RPL24_like_ExpressionFly[RPL24_like_ModuleGenes],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL24_like_ExpressionHuman[RPL24_like_ModuleGenes] and RPL24_like_ExpressionFly[RPL24_like_ModuleGenes]
## S = 137223872, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3482152
```

```r
RPL24_like_shared_module<- data.frame(cbind(RPL24_like_ExpressionHuman = RPL24_like_ExpressionHuman[RPL24_like_ModuleGenes], RPL24_like_ExpressionFly= RPL24_like_ExpressionFly[RPL24_like_ModuleGenes]))

ggplot(RPL24_like_shared_module, aes(RPL24_like_ExpressionHuman, RPL24_like_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-80-1.png" style="display: block; margin: auto;" />

### The _RpL30_- group 
Merge both expression list and metadata data (important to compare the gene co-expression network)

```r
RPL30_treatments<- Fly_treatments %>% dplyr::filter(condition %in% c("RPL30", "control"))
RPL30_counts <- Minute_counts %>% dplyr::select(fly_gene_ID, c(RPL30_treatments$sample))
```

To compare two very distinct data sets, such as humans and flies, it is necessary to analyse the data in the same way. Therefore, we need to combine both Minute and Human data into a vector list.


```r
nSets = 2; #number of Data sets we will work (one for Minute and other for the BioProject PRJNA38870)
setLabels = c("PRJNA497590", "RpL30")

RPL30_gene_expression = vector(mode = "list", length = nSets)

RPL30_gene_expression[[1]] = list(data = as.data.frame(t(PRJNA497590_counts[,-1])));
names(RPL30_gene_expression[[1]]$data) = PRJNA497590_counts$fly_gene_ID;
RPL30_gene_expression[[2]] = list(data = as.data.frame(t(RPL30_counts[,-1])));
names(RPL30_gene_expression[[2]]$data) = RPL30_counts$fly_gene_ID;

##check the list size
exprSize = checkSets(RPL30_gene_expression)
exprSize
```
Removing gene outliers

```r
gsg = goodSamplesGenesMS(RPL30_gene_expression, verbose = 2);
gsg$allOK
if (!gsg$allOK)
{for (set in 1:exprSize$nSets)
{# Remove the offending genes and samples
   RPL30_gene_expression[[set]]$data = RPL30_gene_expression[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
} # Update exprSize
   exprSize = checkSets(RPL30_gene_expression)
   gsg2 = goodSamplesGenesMS(RPL30_gene_expression, verbose = 2); #check again
   gsg2$allOK
}
```
Call the WGCNA funtion that will estimate the soft threshold and the gene co-expression network between humans and treatment groups withing the minute data

```r
source("Function_WGCNA_gene_coexpresison_network.R")
```
Build a gene coexpresison network using WGCNA package and the function build on: "Function_WGCNA_gene_coexpresison_network.R" 
This function will take a while it also runs in silence, to adjust change this alters "verbose = 0" to verbose = 1 or 2

```r
softpower = softpower(RPL30_gene_expression)
RPL30_softpower = softpower

RPL30_human_net<-network(RPL30_gene_expression[[1]]$data)
RPL30_fly_net<- network(RPL30_gene_expression[[2]]$data)
```
#### 4- General expression and module preservation

To visualise the general preservation of the network, estimate the absolute value of the correlation matrix to a a power (soft-thresholding with the power adjacency function)

```r
# Human network:
RPL30_AdjMatHuman =  abs(cor(RPL30_gene_expression[[1]]$data,use="p"))^RPL30_softpower

# Fly network:
RPL30_AdjMatFly =  abs(cor(RPL30_gene_expression[[2]]$data,use="p"))^RPL30_softpower

## Calculation of the whole network connectivity k:
RPL30_ConnectivityHuman=as.vector(apply(RPL30_AdjMatHuman,2,sum, na.rm=T))

RPL30_ConnectivityFly =as.vector(apply(RPL30_AdjMatFly,2,sum, na.rm=T))

## Scaling k to lie between 0 and 1:
RPL30_ConnectivityHuman=RPL30_ConnectivityHuman %>% psycho::standardize() 
RPL30_ConnectivityFly=RPL30_ConnectivityFly %>% psycho::standardize() 
```

Comparing gene expression human and the Minute project

```r
## Comparing gene expression human and fly:
RPL30_ExpressionHuman = apply(RPL30_gene_expression[[1]]$data,2,mean)
RPL30_ExpressionHuman = RPL30_ExpressionHuman/max(RPL30_ExpressionHuman)
RPL30_ExpressionFly = apply(RPL30_gene_expression[[2]]$data,2,mean)
RPL30_ExpressionFly= RPL30_ExpressionFly/max(RPL30_ExpressionFly)
RPL30_cor_data<- data.frame(cbind(RPL30_ExpressionHuman, RPL30_ExpressionFly))

cor.test(RPL30_cor_data$RPL30_ExpressionHuman,RPL30_cor_data$RPL30_ExpressionFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL30_cor_data$RPL30_ExpressionHuman and RPL30_cor_data$RPL30_ExpressionFly
## S = 33565892901, p-value < 0.00000000000000022
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3098882
```

```r
ggplot(RPL30_cor_data, aes(RPL30_ExpressionHuman, RPL30_ExpressionFly)) + 
   geom_point(alpha = 0.3) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-87-1.png" style="display: block; margin: auto;" />

Comparing network connectivity between human and the Minute project

```r
cor.test(RPL30_ConnectivityHuman,RPL30_ConnectivityFly,method="s",use="p") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL30_ConnectivityHuman and RPL30_ConnectivityFly
## S = 48561435557, p-value = 0.8976
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.00158112
```

```r
RPL30_net_data<- data.frame(cbind(RPL30_ConnectivityHuman,RPL30_ConnectivityFly))

ggplot(RPL30_net_data, aes(RPL30_ConnectivityHuman, RPL30_ConnectivityFly)) + 
   geom_point(alpha = 0.3) + # xlim(0, 0.36) + ylim(0, 0.36) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-88-1.png" style="display: block; margin: auto;" />

Identify the genes of interest on module RpL30 in the Minute data

```r
#Isolate RpL30 gene module
RPL30_Module <- RPL30_human_net$colors[which(colnames(RPL30_gene_expression[[1]]$data) == "FBgn0086710")]
RPL30_ModuleGenes <- colnames(RPL30_gene_expression[[1]]$data)[RPL30_human_net$colors == RPL30_Module]
```

####  5- Module preservation 


```r
setLabels <- c("ZikaV", "Experiment") 
ZikaV_RPL30_expression  <- list(Experiment = list(data = RPL30_gene_expression[[2]]$data), ZikaV = list(data = RPL30_gene_expression[[1]]$data))
ZikaV_RPL30_multiColor <- list(Experiment = labels2colors(RPL30_fly_net$colors))
```


```r
#Module preservation
npermutation<- 500
ZikaV_RPL30_mp<- modulePreservation(ZikaV_RPL30_expression , ZikaV_RPL30_multiColor, nPermutations = npermutation,
                                          referenceNetworks =1,randomSeed = 1, quickCor = 0, verbose = 3)
#save(ZikaV_RPL30_mp, file = "RDS/ZikaV_RPL30_mp.RData");

#_P-value_
ZikaV_RPL30_p.value=exp(ZikaV_RPL30_mp$preservation$log.p$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL30_Module),"log.psummary.pres"])
ZikaV_RPL30_Zsummary = ZikaV_RPL30_mp$preservation$Z$ref.Experiment$inColumnsAlsoPresentIn.ZikaV[labels2colors(RPL30_Module),"Zsummary.pres"]

ZikaV_mp= rbind(ZikaV_mp, "RpL30" = c(ZikaV_RPL30_Zsummary,ZikaV_RPL30_p.value))
```
Comparing correlation of the RpL30 module between humans and the Minute project


```r
cor.test(RPL30_ExpressionHuman[RPL30_ModuleGenes],RPL30_ExpressionFly[RPL30_ModuleGenes ],method="s",use="p")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  RPL30_ExpressionHuman[RPL30_ModuleGenes] and RPL30_ExpressionFly[RPL30_ModuleGenes]
## S = 3242029, p-value = 0.000000000000000719
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.4280773
```

```r
RPL30_shared_module<- data.frame(cbind(RPL30_ExpressionHuman = RPL30_ExpressionHuman[RPL30_ModuleGenes], RPL30_ExpressionFly= RPL30_ExpressionFly[RPL30_ModuleGenes]))

ggplot(RPL30_shared_module, aes(RPL30_ExpressionHuman, RPL30_ExpressionFly)) + 
   geom_point(alpha = 0.3) +  xlim(0, 0.12) + ylim(0, 0.25) + 
   geom_smooth(aes(group=1), method = "lm", se=F, color = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<img src="Project_PRJNA497590/figures/unnamed-chunk-92-1.png" style="display: block; margin: auto;" />

## Module preservation


```r
library(ztable)
```

```
## Welcome to package ztable ver 0.2.0
```

```r
options(ztable.type="html")
ztable(as.data.frame(ZikaV_mp),digits = 7) %>% makeHeatmap(mycolor=gradientColor(high="blue",mid="white",low="red"),cols=2, reverse = F) %>%  makeHeatmap(palette="Blues",cols=1) %>% print(caption="Zika Virus and Minute module preservation analysis")
```

<head><style>
        table {
              font-family: times ;
color:  black ;
text-align: right;}
        th {
              padding: 1px 1px 5px 5px;
	        }
        td {
             padding: 1px 1px 5px 5px; }
      </style></head><table align="center" style="border-collapse: collapse; caption-side:top; font-size:11pt;"><caption style="text-align:center;">Zika Virus and Minute module preservation analysis</caption><tr>
<th style="border-left: 0px solid black;background-color: #FFFFFF;border-top: 2px solid gray;border-bottom: 1px solid gray;">&nbsp;</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">Zsummary</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">p.value</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">RpL3</td>
<td align="right" style="border-left: 0px solid black;background-color: #9ECAE1;">0.0940722</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;background-color: #B38BFF;color: #FFFFFF;">0.7100884</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpL14</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #F7FBFF;">-1.4252926</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #0000FF;color: #FFFFFF;">0.9288878</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpS13</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #08306B;color: #FFFFFF;">2.1782467</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FF4124;">0.1488207</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpL19</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #6BAED6;color: #FFFFFF;">0.2974175</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #E3D0FF;color: #FFFFFF;">0.6009009</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpS19b</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #08519C;color: #FFFFFF;">1.8472781</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FF6140;">0.2155471</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpL24-like</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #08306B;color: #FFFFFF;">2.2415883</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FF0000;">0.1154291</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">RpL30</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #DEEBF7;">-0.8212371</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #7145FF;color: #FFFFFF;">0.8375188</td>
</tr>
<tr>
<td colspan="3" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>


## Session info

```r
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 3.6.1 (2019-07-05)
##  os       macOS Catalina 10.15.5      
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_AU.UTF-8                 
##  ctype    en_AU.UTF-8                 
##  tz       Australia/Melbourne         
##  date     2020-07-05                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package         * version  date       lib source                             
##  acepack           1.4.1    2016-10-29 [1] CRAN (R 3.6.0)                     
##  AnnotationDbi     1.48.0   2019-10-29 [1] Bioconductor                       
##  assertthat        0.2.1    2019-03-21 [2] CRAN (R 3.6.0)                     
##  backports         1.1.8    2020-06-17 [1] CRAN (R 3.6.2)                     
##  base64enc         0.1-3    2015-07-28 [1] CRAN (R 3.6.0)                     
##  Biobase           2.46.0   2019-10-29 [1] Bioconductor                       
##  BiocGenerics    * 0.32.0   2019-10-29 [1] Bioconductor                       
##  BiocManager       1.30.10  2019-11-16 [1] CRAN (R 3.6.0)                     
##  BiocParallel      1.20.1   2019-12-21 [1] Bioconductor                       
##  bit               1.1-15.2 2020-02-10 [1] CRAN (R 3.6.0)                     
##  bit64             0.9-7    2017-05-08 [1] CRAN (R 3.6.0)                     
##  blob              1.2.1    2020-01-20 [1] CRAN (R 3.6.0)                     
##  broom             0.5.6    2020-04-20 [1] CRAN (R 3.6.2)                     
##  callr             3.4.3    2020-03-28 [1] CRAN (R 3.6.2)                     
##  cellranger        1.1.0    2016-07-27 [1] CRAN (R 3.6.0)                     
##  checkmate         2.0.0    2020-02-06 [1] CRAN (R 3.6.0)                     
##  chorddiag       * 0.1.2    2020-03-05 [1] Github (mattflor/chorddiag@86b1652)
##  circlize        * 0.4.10   2020-06-15 [1] CRAN (R 3.6.1)                     
##  cli               2.0.2    2020-02-28 [2] CRAN (R 3.6.0)                     
##  cluster           2.1.0    2019-06-19 [2] CRAN (R 3.6.1)                     
##  clusterProfiler * 3.14.3   2020-01-08 [1] Bioconductor                       
##  codetools         0.2-16   2018-12-24 [2] CRAN (R 3.6.1)                     
##  colorspace        1.4-1    2019-03-18 [2] CRAN (R 3.6.0)                     
##  cowplot           1.0.0    2019-07-11 [1] CRAN (R 3.6.0)                     
##  crayon            1.3.4    2017-09-16 [2] CRAN (R 3.6.0)                     
##  data.table        1.12.8   2019-12-09 [1] CRAN (R 3.6.0)                     
##  DBI               1.1.0    2019-12-15 [1] CRAN (R 3.6.0)                     
##  dbplyr            1.4.4    2020-05-27 [1] CRAN (R 3.6.2)                     
##  desc              1.2.0    2018-05-01 [2] CRAN (R 3.6.0)                     
##  devtools          2.3.0    2020-04-10 [1] CRAN (R 3.6.2)                     
##  digest            0.6.25   2020-02-23 [2] CRAN (R 3.6.0)                     
##  DO.db             2.9      2020-04-12 [1] Bioconductor                       
##  doParallel        1.0.15   2019-08-02 [1] CRAN (R 3.6.0)                     
##  DOSE              3.12.0   2019-10-29 [1] Bioconductor                       
##  dplyr           * 1.0.0    2020-05-29 [1] CRAN (R 3.6.2)                     
##  dynamicTreeCut  * 1.63-1   2016-03-11 [1] CRAN (R 3.6.0)                     
##  ellipsis          0.3.1    2020-05-15 [2] CRAN (R 3.6.2)                     
##  enrichplot        1.6.1    2019-12-16 [1] Bioconductor                       
##  europepmc         0.4      2020-05-31 [1] CRAN (R 3.6.1)                     
##  evaluate          0.14     2019-05-28 [1] CRAN (R 3.6.0)                     
##  fansi             0.4.1    2020-01-08 [2] CRAN (R 3.6.0)                     
##  farver            2.0.3    2020-01-16 [2] CRAN (R 3.6.0)                     
##  fastcluster     * 1.1.25   2018-06-07 [1] CRAN (R 3.6.0)                     
##  fastmatch         1.1-0    2017-01-28 [1] CRAN (R 3.6.0)                     
##  fgsea             1.12.0   2019-10-29 [1] Bioconductor                       
##  flashClust      * 1.01-2   2012-08-21 [1] CRAN (R 3.6.0)                     
##  flextable         0.5.10   2020-05-15 [1] CRAN (R 3.6.2)                     
##  forcats         * 0.5.0    2020-03-01 [1] CRAN (R 3.6.0)                     
##  foreach           1.5.0    2020-03-30 [1] CRAN (R 3.6.2)                     
##  foreign           0.8-76   2020-03-03 [2] CRAN (R 3.6.0)                     
##  Formula           1.2-3    2018-05-03 [1] CRAN (R 3.6.0)                     
##  fs                1.4.1    2020-04-04 [1] CRAN (R 3.6.1)                     
##  gdtools           0.2.2    2020-04-03 [1] CRAN (R 3.6.2)                     
##  generics          0.0.2    2018-11-29 [1] CRAN (R 3.6.0)                     
##  ggforce           0.3.1    2019-08-20 [1] CRAN (R 3.6.0)                     
##  ggplot2         * 3.3.2    2020-06-19 [2] CRAN (R 3.6.2)                     
##  ggplotify         0.0.5    2020-03-12 [1] CRAN (R 3.6.0)                     
##  ggraph          * 2.0.3    2020-05-20 [1] CRAN (R 3.6.2)                     
##  ggrepel           0.8.2    2020-03-08 [1] CRAN (R 3.6.0)                     
##  ggridges          0.5.2    2020-01-12 [1] CRAN (R 3.6.0)                     
##  GlobalOptions     0.1.2    2020-06-10 [1] CRAN (R 3.6.2)                     
##  glue              1.4.1    2020-05-13 [1] CRAN (R 3.6.2)                     
##  GO.db             3.10.0   2020-03-06 [1] Bioconductor                       
##  GOSemSim          2.12.1   2020-03-19 [1] Bioconductor                       
##  graphlayouts      0.7.0    2020-04-25 [1] CRAN (R 3.6.2)                     
##  gridExtra         2.3      2017-09-09 [1] CRAN (R 3.6.0)                     
##  gridGraphics      0.5-0    2020-02-25 [1] CRAN (R 3.6.0)                     
##  gtable            0.3.0    2019-03-25 [2] CRAN (R 3.6.0)                     
##  haven             2.3.1    2020-06-01 [1] CRAN (R 3.6.2)                     
##  Hmisc             4.4-0    2020-03-23 [1] CRAN (R 3.6.0)                     
##  hms               0.5.3    2020-01-08 [1] CRAN (R 3.6.0)                     
##  htmlTable         2.0.0    2020-06-21 [1] CRAN (R 3.6.1)                     
##  htmltools         0.5.0    2020-06-16 [1] CRAN (R 3.6.2)                     
##  htmlwidgets       1.5.1    2019-10-08 [1] CRAN (R 3.6.0)                     
##  httr              1.4.1    2019-08-05 [1] CRAN (R 3.6.0)                     
##  igraph          * 1.2.5    2020-03-19 [1] CRAN (R 3.6.0)                     
##  impute            1.60.0   2019-10-29 [1] Bioconductor                       
##  IRanges         * 2.20.2   2020-01-13 [1] Bioconductor                       
##  iterators         1.0.12   2019-07-26 [1] CRAN (R 3.6.0)                     
##  jpeg              0.1-8.1  2019-10-24 [1] CRAN (R 3.6.0)                     
##  jsonlite          1.6.1    2020-02-02 [1] CRAN (R 3.6.0)                     
##  kableExtra      * 1.1.0    2019-03-16 [1] CRAN (R 3.6.0)                     
##  knitr           * 1.28     2020-02-06 [1] CRAN (R 3.6.0)                     
##  labeling          0.3      2014-08-23 [2] CRAN (R 3.6.0)                     
##  lattice           0.20-41  2020-04-02 [2] CRAN (R 3.6.1)                     
##  latticeExtra      0.6-29   2019-12-19 [1] CRAN (R 3.6.0)                     
##  lifecycle         0.2.0    2020-03-06 [2] CRAN (R 3.6.0)                     
##  lubridate         1.7.9    2020-06-08 [1] CRAN (R 3.6.2)                     
##  magrittr          1.5      2014-11-22 [2] CRAN (R 3.6.0)                     
##  MASS              7.3-51.6 2020-04-26 [2] CRAN (R 3.6.2)                     
##  Matrix            1.2-18   2019-11-27 [2] CRAN (R 3.6.0)                     
##  matrixStats       0.56.0   2020-03-13 [1] CRAN (R 3.6.0)                     
##  memoise           1.1.0    2017-04-21 [1] CRAN (R 3.6.0)                     
##  mgcv              1.8-31   2019-11-09 [2] CRAN (R 3.6.0)                     
##  modelr            0.1.8    2020-05-19 [1] CRAN (R 3.6.2)                     
##  munsell           0.5.0    2018-06-12 [2] CRAN (R 3.6.0)                     
##  nlme              3.1-148  2020-05-24 [2] CRAN (R 3.6.2)                     
##  nnet              7.3-14   2020-04-26 [2] CRAN (R 3.6.2)                     
##  officer           0.3.11   2020-05-18 [1] CRAN (R 3.6.2)                     
##  pillar            1.4.4    2020-05-05 [2] CRAN (R 3.6.2)                     
##  pkgbuild          1.0.8    2020-05-07 [2] CRAN (R 3.6.2)                     
##  pkgconfig         2.0.3    2019-09-22 [2] CRAN (R 3.6.0)                     
##  pkgload           1.1.0    2020-05-29 [2] CRAN (R 3.6.2)                     
##  plyr              1.8.6    2020-03-03 [2] CRAN (R 3.6.0)                     
##  png               0.1-7    2013-12-03 [1] CRAN (R 3.6.0)                     
##  polyclip          1.10-0   2019-03-14 [1] CRAN (R 3.6.0)                     
##  preprocessCore    1.48.0   2019-10-29 [1] Bioconductor                       
##  prettyunits       1.1.1    2020-01-24 [1] CRAN (R 3.6.0)                     
##  processx          3.4.2    2020-02-09 [1] CRAN (R 3.6.0)                     
##  progress          1.2.2    2019-05-16 [1] CRAN (R 3.6.0)                     
##  ps                1.3.3    2020-05-08 [1] CRAN (R 3.6.2)                     
##  psycho          * 0.5.0    2020-01-22 [1] CRAN (R 3.6.0)                     
##  purrr           * 0.3.4    2020-04-17 [1] CRAN (R 3.6.2)                     
##  qvalue            2.18.0   2019-10-29 [1] Bioconductor                       
##  R6                2.4.1    2019-11-12 [2] CRAN (R 3.6.0)                     
##  RColorBrewer      1.1-2    2014-12-07 [1] CRAN (R 3.6.0)                     
##  Rcpp              1.0.4.6  2020-04-09 [1] CRAN (R 3.6.1)                     
##  readr           * 1.3.1    2018-12-21 [1] CRAN (R 3.6.0)                     
##  readxl            1.3.1    2019-03-13 [1] CRAN (R 3.6.0)                     
##  remotes           2.1.1    2020-02-15 [1] CRAN (R 3.6.0)                     
##  reprex            0.3.0    2019-05-16 [1] CRAN (R 3.6.0)                     
##  reshape2          1.4.4    2020-04-09 [2] CRAN (R 3.6.1)                     
##  rlang             0.4.6    2020-05-02 [2] CRAN (R 3.6.2)                     
##  rmarkdown         2.3      2020-06-18 [1] CRAN (R 3.6.1)                     
##  rpart             4.1-15   2019-04-12 [2] CRAN (R 3.6.1)                     
##  rprojroot         1.3-2    2018-01-03 [2] CRAN (R 3.6.0)                     
##  RSQLite           2.2.0    2020-01-07 [1] CRAN (R 3.6.0)                     
##  rstudioapi        0.11     2020-02-07 [1] CRAN (R 3.6.0)                     
##  rvcheck           0.1.8    2020-03-01 [1] CRAN (R 3.6.0)                     
##  rvest             0.3.5    2019-11-08 [1] CRAN (R 3.6.0)                     
##  S4Vectors       * 0.24.4   2020-04-09 [1] Bioconductor                       
##  scales            1.1.1    2020-05-11 [2] CRAN (R 3.6.2)                     
##  sessioninfo       1.1.1    2018-11-05 [1] CRAN (R 3.6.0)                     
##  shape             1.4.4    2018-02-07 [1] CRAN (R 3.6.0)                     
##  stringi           1.4.6    2020-02-17 [1] CRAN (R 3.6.0)                     
##  stringr         * 1.4.0    2019-02-10 [1] CRAN (R 3.6.0)                     
##  survival          3.2-3    2020-06-13 [1] CRAN (R 3.6.2)                     
##  systemfonts       0.2.3    2020-06-09 [1] CRAN (R 3.6.2)                     
##  testthat          2.3.2    2020-03-02 [2] CRAN (R 3.6.0)                     
##  tibble          * 3.0.1    2020-04-20 [1] CRAN (R 3.6.2)                     
##  tidygraph         1.2.0    2020-05-12 [1] CRAN (R 3.6.2)                     
##  tidyr           * 1.1.0    2020-05-20 [1] CRAN (R 3.6.2)                     
##  tidyselect        1.1.0    2020-05-11 [1] CRAN (R 3.6.2)                     
##  tidyverse       * 1.3.0    2019-11-21 [1] CRAN (R 3.6.0)                     
##  triebeard         0.3.0    2016-08-04 [1] CRAN (R 3.6.0)                     
##  tweenr            1.0.1    2018-12-14 [1] CRAN (R 3.6.0)                     
##  tximport        * 1.14.2   2020-03-22 [1] Bioconductor                       
##  urltools          1.7.3    2019-04-14 [1] CRAN (R 3.6.0)                     
##  usethis           1.6.1    2020-04-29 [1] CRAN (R 3.6.2)                     
##  uuid              0.1-4    2020-02-26 [1] CRAN (R 3.6.0)                     
##  vctrs             0.3.1    2020-06-05 [1] CRAN (R 3.6.2)                     
##  viridis         * 0.5.1    2018-03-29 [1] CRAN (R 3.6.0)                     
##  viridisLite     * 0.3.0    2018-02-01 [2] CRAN (R 3.6.0)                     
##  webshot           0.5.2    2019-11-22 [1] CRAN (R 3.6.0)                     
##  WGCNA           * 1.69     2020-02-28 [1] CRAN (R 3.6.1)                     
##  withr             2.2.0    2020-04-20 [1] CRAN (R 3.6.2)                     
##  xfun              0.14     2020-05-20 [1] CRAN (R 3.6.2)                     
##  xml2              1.3.2    2020-04-23 [1] CRAN (R 3.6.2)                     
##  yaml              2.2.1    2020-02-01 [1] CRAN (R 3.6.0)                     
##  zip               2.0.4    2019-09-01 [1] CRAN (R 3.6.0)                     
##  ztable          * 0.2.0    2018-07-22 [1] CRAN (R 3.6.0)                     
## 
## [1] /Users/mvel0008/Library/R/3.6/library
## [2] /Library/Frameworks/R.framework/Versions/3.6/Resources/library
```
