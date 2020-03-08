kwithin_dat <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  kWithin <- sapply(unlist(tgenes), function(x) genes$kWithin[match(x, genes$ID)])
  if(class(kWithin) == 'factor'){
    kWithin <- gsub(",", ".", gsub("\\.", "", kWithin))
    kWithin <- as.numeric(kWithin)
  }
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), '/')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  kDiff <- sapply(unlist(tgenes), function(x) genes$kDiff[match(x, genes$ID)])
  if(class(kDiff) == 'factor'){
    kDiff <- gsub(",", ".", gsub("\\.", "", kDiff))
    kDiff <- as.numeric(kDiff)
  }
if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), '/')
count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
gene_id<- sapply(unlist(tgenes), function(x) genes$gene_id[match(x, genes$ID)])
if(class(gene_id) == 'factor'){
  gene_id<- gsub(",", ".", gsub("\\.", "", gene_id))
  gene_id<- as.factor(gene_id)
}
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(kWithin[s:e], function(x) ifelse(x > 0, 1, -1))
    zsc <- c(zsc, sum(value) / sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(genes = as.character(unlist(tgenes)),category = rep(as.character(terms$category), count),gene_id = gene_id, 
                     term = rep(as.character(terms$term), count),count = rep(count, count),  kDiff  = kDiff, kWithin = kWithin,
                     adj_pval = rep(terms$adj_pval, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(genes = as.character(unlist(tgenes)),category = rep(as.character(terms$category), count),gene_id = gene_id, 
                     term = rep(as.character(terms$term), count),count = rep(count, count),  kDiff  = kDiff, kWithin = kWithin,
                     adj_pval = rep(terms$adj_pval, count), stringsAsFactors = FALSE)
  }
  return(df)
}



dat_chord <- function(data, gene, process){
  id <- term <- kWithin <- BPprocess <- NULL
  
  colnames(data) <- tolower(colnames(data))
  if (missing(genes)){
    if (is.null(data$kWithin)){
      genes <- as.character(unique(data$genes))
    }else{
      genes <- subset(data, !duplicated(genes), c(genes, kWithin))
    }
  }else{
    if(is.vector(genes)){
      genes <- as.character(genes) 
    }else{
      if(class(genes[, 2]) != 'numeric') genes[, 2] <- as.numeric(levels(genes[, 2]))[genes[, 2]]
      genes[, 1] <- as.character(genes[, 1])
      colnames(genes) <- c('genes', 'kWithin')
    }
  }
  if (missing(process)){
    process <- as.character(unique(data$term))
  }else{
    if(class(process) != 'character') process <- as.character(process)
  }
  if (strsplit(process[1],':')[[1]][1] == 'GO'){
    subData <- subset(data, id%in%process)
    colnames(subData)[which(colnames(subData) == 'id')] <- 'BPprocess'
  }else{
    subData <- subset(data, term%in%process)
    colnames(subData)[which(colnames(subData) == 'term')] <- 'BPprocess'
  }
  
  if(is.vector(genes)){
    M <- genes[genes%in%unique(subData$genes)]
    mat <- matrix(0, ncol = length(process), nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- process
    for (p in 1:length(process)){
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
    }
  }else{
    genes <- subset(genes, genes %in% unique(subData$genes))
    N <- length(process) + 1
    M <- genes[,1] 
    mat <- matrix(0, ncol = N, nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- c(process, 'kWithin') 
    mat[,N] <- genes[,2]
    for (p in 1:(N-1)){
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
    }
  }
  return(mat)
}

