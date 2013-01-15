#
# Contains functions for performing PathNet Contextual Analysis
#

# do.contextual.analysis
#   Main function for performing contextual analysis
# Arguments
#   DirectEvidence - direct evidence data
#   Adjacency - Adjacency matrix
#   gene_ID - gene IDs
#   pathway - pathway data
#   pathwayNames - name of pathways
#   n_perm - number of permutations to run
#   use_sig_pathways - whether to use significant pathway analysis
#   ea.results - enrichment analysis results
# Returns contextual association between pathways and pathway overlaps
#   list(conn_p_value, overlap)
#     conn_p_value - matrix of contextual p values
#     overlap - pathway overlap matrix
#
do.contextual.analysis <- function(DirectEvidence, Adjacency, gene_ID, pathway, pathwayNames, n_perm, use_sig_pathways, ea.results=NULL) {
  # Setup initial variables
  Adjacency <- Adjacency[rownames(Adjacency) %in% gene_ID,rownames(Adjacency) %in% gene_ID]  
  observe <- rep(FALSE, nrow(Adjacency))
  names(observe) <- rownames(Adjacency)
  diag(Adjacency) <- 0
  path_genes_all  <- unique(c(pathway[,1], pathway[,2]))
  
  # Setup the pathway names we want to use
  pathwayNames_significant <- construct.significant.pathways(use_sig_pathways, pathwayNames, ea.results)
  message("Number of significant pathways to use for contextual analysis:", length(pathwayNames_significant))
  
  Conn_p_value <- array(1, dim <- c(length(pathwayNames_significant), length(pathwayNames_significant)))
  overlap      <- array(1, dim <- c(length(pathwayNames_significant), length(pathwayNames_significant)))

  # If we have pathways we can process
  if (length(pathwayNames_significant) > 1)
  {
    #
    # Main contextual processing loop
    #
    for (i in 1:length(pathwayNames_significant))
    {
      message(paste("Performing contextual analysis of",pathwayNames_significant[i],"pathway"))
      
      # Determine common path for current pathway
      common_path_results <- determine.common.path(pathwayNames_significant[i],
                                                   pathway, DirectEvidence, 
                                                   Adjacency, gene_ID)
      path_genes1                 <- common_path_results$path_genes1
      path_genes_common1          <- common_path_results$path_genes_common1
      DirectEvidence_path_common1 <- common_path_results$DirectEvidence_path_common1
      
      # Mail loop to calculate contextual intersection and overlaps
      for (j in 1:length(pathwayNames_significant))
      {
        ind2 <- pathway[,3] == pathwayNames_significant[j]
        path_genes2 <- unique(c(pathway[ind2,1],pathway[ind2,2]))
        path_genes2 <- path_genes2[path_genes2 != "0"]
        path_genes_common2 <- path_genes2 
        
        # find intersections
        if (is.null(path_genes_common2) == FALSE)
        {
          path_genes_common2 <- intersect(path_genes_common2, as.matrix(rownames(Adjacency)))
          if (is.null(path_genes_common2) == FALSE)
          {
            path_genes_common2 <- intersect(path_genes_common2,gene_ID)
          }
        }
        
        # Determine whether to use this path
        use.path.results <- determine.use.path(path_genes_common1, path_genes_common2, Adjacency)
        Adjacency_path   <- use.path.results$Adjacency_path
        
        if (use.path.results$use_this_path == 1)
        {
          degree_Adjacency <- rowSums(Adjacency_path)
          DirectEvidence_path_common2 <- sapply(path_genes_common2, function(p) { DirectEvidence[gene_ID %in% p] })
          
          Conn <- DirectEvidence_path_common1 %*% Adjacency_path 
          Conn <- Conn %*% DirectEvidence_path_common2
          Conn_rand <- sapply(1:n_perm, function(i) {
            DirectEvidence_path_rand1 <- sample(DirectEvidence, length(DirectEvidence_path_common1), replace = FALSE)
            DirectEvidence_path_rand2 <- sample(DirectEvidence, length(DirectEvidence_path_common2), replace = FALSE)
            temp <- DirectEvidence_path_rand1 %*% Adjacency_path
            temp <- temp %*% DirectEvidence_path_rand2
            temp
          })
          temp <- rep(0,n_perm)
          Conn <- rep(Conn,n_perm)
          temp[Conn_rand >= Conn] <- 1
          Conn_p_value[i,j] <- sum(temp) / n_perm
        }
        
        # Construct overlap
        overlap[i,j] <- phyper(q = length(intersect(path_genes1,path_genes2))-1,
                               m = length(path_genes1),
                               n = length(path_genes_all) - length(path_genes1),
                               k = length(path_genes2), lower.tail = FALSE)		
      }
    }
  }
  else {
    warning("Number of significance pathways is less than two")
  }
  
  #
  # Construct results and return contextual association 
  # between pathways and pathway overlaps
  #
  diag(Conn_p_value) <- 0
  diag(overlap) <- 0
  rownames(Conn_p_value) <- pathwayNames_significant 
  colnames(Conn_p_value) <- pathwayNames_significant 
  rownames(overlap) <- pathwayNames_significant 
  colnames(overlap) <- pathwayNames_significant
  list(conn_p_value=Conn_p_value, overlap=overlap)
}

# Determines which pathways to use for contextual analysis
# Arguments
#   use_sig_pathways - whether to only use significant pathways
#   pathwayNames - names of pathways
#   ea.results - Enrichment Analysis results containing correct p values for significance
#
construct.significant.pathways <- function(use_sig_pathways, pathwayNames, ea.results) {
  if (use_sig_pathways == TRUE)
  {
    rankwise <- ea.results$rankwise
    path_pVal_p.Comb_corrected_Bon <- ea.results$result$p_PathNet_FWER
    
    pathwayNames_rankwise <- pathwayNames[rankwise]
    ind <- path_pVal_p.Comb_corrected_Bon[rankwise] <= 0.05
    pathwayNames_significant <- pathwayNames_rankwise[ind]
  }
  else
  {
    # Just sort pathways
    pathwayNames_significant <- sort(pathwayNames)
  }
  pathwayNames_significant
}

# Determines common pathway during analysis
determine.common.path <- function(pathwayName, pathway, DirectEvidence, Adjacency, gene_ID) {
  ind1 <- pathway[,3] == pathwayName
  path_genes1 <- unique(c(pathway[ind1,1],pathway[ind1,2]))
  path_genes1 <- path_genes1[path_genes1 != "0"]
  path_genes_common1 <- path_genes1
  DirectEvidence_path_common1 = NULL
  
  if (is.null(path_genes_common1) == FALSE)
  {
    path_genes_common1 <- intersect(path_genes_common1, as.matrix(rownames(Adjacency)))
    if (is.null(path_genes_common1) == FALSE)
    {
      path_genes_common1 <- intersect(path_genes_common1, gene_ID)
      if (length(path_genes_common1) > 1)
      {
        DirectEvidence_path_common1 <- sapply(path_genes_common1, function(p) { DirectEvidence[gene_ID %in% p] })
      }
    }
  }
  list(path_genes1=path_genes1, 
       path_genes_common1=path_genes_common1, 
       DirectEvidence_path_common1=DirectEvidence_path_common1)
}

# Determines whether to use the specified path, constructs correspoding adjacency path
determine.use.path <- function(path_genes_common1, path_genes_common2, Adjacency) {
  use_this_path <- 0
  Adjacency_path <- NULL
  if (is.null(path_genes_common1) == FALSE && is.null(path_genes_common2) == FALSE && 
        length(path_genes_common1) > 1 && length(path_genes_common2) > 1) {
    Adjacency_path <- Adjacency[as.matrix(rownames(Adjacency)) %in% path_genes_common1,
                                as.matrix(rownames(Adjacency)) %in% path_genes_common2]  
    use_this_path <- 1	
  }
  list(use_this_path=use_this_path, Adjacency_path=Adjacency_path)
}
