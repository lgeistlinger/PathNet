# Main function for performing enrichment analysis
# Returns
#   list(
#     result - main enrichment results
#     combined.evidence - combined evidence (indirect and direct)
#     rankwise rankwise ordering)
do.enrichment.analysis <- function(DirectEvidence, Adjacency, observe, gene_ID, 
                                   pathway, pathwayNames, n_perm, threshold) {
  message("Starting enrichment analysis")
  
  # Create variables for direct evidence
  DirectEvidence_rowwise <- sapply(1:nrow(Adjacency),
                                   function(i) { 
                                     DirectEvidence[gene_ID %in% names(observe[i])]})
  
  pDirectEvidence <- 10^(-1*DirectEvidence)
  pComb <- pDirectEvidence
  
  # Create variables for indirect evidence
  IndirectEvidence_p_value <- rep(NA, nrow(Adjacency))
  IndirectEvidence_return <- rep(NA,length(gene_ID))  # IndirectEvidence_return is calculated for all genes, including the ones not connected
  IndirectEvidence_rand <- rep(NA,nrow(Adjacency)*n_perm)
  IndirectEvidence_rand <- matrix(IndirectEvidence_rand, nrow = nrow(Adjacency), ncol = n_perm)
  degree_Adjacency <- rowSums(Adjacency)
  n_path_genes <- rep(NA,length(pathwayNames))
  IndirectEvidence_return <- rep(NA,length(gene_ID))
  
  # Construct pathway links
  pathway_links <- construct.pathway_links(pathway, pathwayNames)
  path_genes_all  <- unique(c(pathway[,1], pathway[,2]))
  path_genes_all  <- path_genes_all[path_genes_all != 0]
  
  # perform random sampling for indirect evidence
  IndirectEvidence <- Adjacency %*% DirectEvidence_rowwise
  IndirectEvidence_rand <- replicate(n_perm, 
                                     Adjacency %*% sample(DirectEvidence, nrow(Adjacency), replace=FALSE),
                                     simplify="matrix")
  
  
  #
  # Calculate indirect evidence p values and combined probabilities
  #
  for ( i in 1:nrow(Adjacency))
  {
    if (degree_Adjacency[i] > 0)   # For isolated genes in KEGG pathway only expression p-value will be used
    {
      temp <- rep(0,n_perm)
      if (IndirectEvidence[i] > 0)
        temp[(IndirectEvidence_rand[i,]>=IndirectEvidence[i])] <- 1
      IndirectEvidence_p_value[i] <- sum(temp)/n_perm
      if (IndirectEvidence_p_value[i] >1)
        IndirectEvidence_p_value[i] = 1			
      if (IndirectEvidence_p_value[i] < 1/n_perm )
        IndirectEvidence_p_value[i] = 1/n_perm 			
      IndirectEvidence_return[gene_ID %in% names(observe[i])] <- IndirectEvidence_p_value[i]  #   It has the same size as original data
      product <- IndirectEvidence_p_value[i]*(10^(-1*DirectEvidence_rowwise[i]))
      pComb[gene_ID %in% names(observe[i])] <- product - product * log(product)
      rm(temp,product)
    }
  }
  
  combined.evidence <- data.frame(gene_ID = gene_ID,
                                  pDirectEvidence = pDirectEvidence,
                                  pIndirectEvidence = IndirectEvidence_return,
                                  pCombinedEvidence = pComb)  
  
  # ****************************************Calculation of pathway enrichment***********************************************
  path_genes_all <- intersect(path_genes_all,gene_ID)
  
  sig_genes_p.Comb <- intersect(gene_ID[pComb < threshold],path_genes_all)
  sig_genes_p.Original <- intersect(gene_ID[pDirectEvidence < threshold],path_genes_all)
  
  pathway_sig_p.Comb <- rep(NA,length(pathwayNames))
  pathway_sig_p.Original <- rep(NA,length(pathwayNames))
  path_pVal_p.Comb <- rep(NA,length(pathwayNames))
  path_pVal_p.Original <- rep(NA,length(pathwayNames))
  
  #
  # Main pathway calculation loop
  #
  for (i in 1:length(pathwayNames))
  {
    ind <- pathway[,3] == pathwayNames[i]
    path_genes <- intersect(c(pathway[ind,1],pathway[ind,2]),path_genes_all)
    path_genes <- path_genes[path_genes != "0"]
    n_path_genes[i] <- length(path_genes)
    pathway_sig_p.Comb[i] <- length(intersect(path_genes,sig_genes_p.Comb))
    pathway_sig_p.Original[i] <- length(intersect(path_genes,sig_genes_p.Original))
    if (length(sig_genes_p.Comb) <= 1)
      path_pVal_p.Comb[i] <- 1
    else {
      path_pVal_p.Comb[i] <- phyper(q = pathway_sig_p.Comb[i]-1, 
                                    m = length(path_genes), 
                                    n = length(path_genes_all) - length(path_genes), 
                                    k = length(sig_genes_p.Comb), lower.tail = FALSE)	
    }
    if (length(sig_genes_p.Original) <= 1) {
      path_pVal_p.Original[i] <- 1
    }
    else {
      path_pVal_p.Original[i] <- phyper(q = pathway_sig_p.Original[i]-1, 
                                        m = length(path_genes), 
                                        n = length(path_genes_all) - length(path_genes), 
                                        k = length(sig_genes_p.Original), lower.tail = FALSE)	
    } 
    message(paste("Completed enrichment analysis of",pathwayNames[i],"pathway"))
  } # end main pathway loop
  
  path_pVal_p.Original_corrected_Bon <- p.adjust(path_pVal_p.Original, method  = "bonferroni")
  path_pVal_p.Comb_corrected_Bon <- p.adjust(path_pVal_p.Comb, method = "bonferroni")
  path_pVal_p.Original_corrected_FDR <- p.adjust(path_pVal_p.Original, method  = "fdr")
  path_pVal_p.Comb_corrected_FDR <- p.adjust(path_pVal_p.Comb, method = "fdr")
  
  
  # ***************************************Constructing the results ***************************************************
  rankwise <- order(path_pVal_p.Comb)
  result <- data.frame(Name = substr(pathwayNames[rankwise],1,35),
                       No_of_Genes    = n_path_genes[rankwise],
                       Sig_Direct     = pathway_sig_p.Original[rankwise],
                       Sig_Combi      = pathway_sig_p.Comb[rankwise],
                       p_Hyper        = path_pVal_p.Original[rankwise],
                       p_Hyper_FWER   = path_pVal_p.Original_corrected_Bon[rankwise],
                       p_PathNet      = path_pVal_p.Comb[rankwise], 
                       p_PathNet_FWER = path_pVal_p.Comb_corrected_Bon[rankwise])
  
  message("Number of significant genes from direct evidence")
  message(length(sig_genes_p.Original))
  message("Number of significant genes from combined evidence")  
  message(length(sig_genes_p.Comb))
  # Return results
  list(result=result, combined.evidence=combined.evidence, rankwise=rankwise)
}


# Identifies number of pathway links for each pathway
construct.pathway_links <- function(pathway, pathwayNames) {
  message(paste("Processing ",length(pathwayNames),"pathways"))
  pathway_links <- sapply(pathwayNames, function(name) {
    idx <- pathway[,1] != 0 & pathway[,2] !=0 & pathway[,3] == name
    length(idx[idx==TRUE])
  })
}