#
# Main PathNet Function
#
# See man page for full description.
# 
# Returns a list structure containing
#   enrichment_results
#   enrichment_combined_evidence
#   conn_p_value
#   pathway_overlap
#   

PathNet <- function(Enrichment_Analysis = TRUE, Contextual_Analysis = FALSE, 
                    DirectEvidence_info, Column_DirectEvidence = 7, 
                    Adjacency, pathway = pathway, n_perm = 2000,
                    threshold = 0.05, use_sig_pathways = FALSE) 
{
  # Created by Bhaskar Dutta @BHSAI
  #
  #    This program is free software: you can redistribute it and/or modify
  #    it under the terms of the GNU General Public License as published by
  #    the Free Software Foundation, either version 3 of the License, or
  #    (at your option) any later version.
  #    
  #    This program is distributed in the hope that it will be useful,
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #    GNU General Public License for more details.
  #    
  #    You should have received a copy of the GNU General Public License
  #    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  #    
  
  # Verify arguments
  if (use_sig_pathways && !Enrichment_Analysis) {
    stop('Using significant pathways requires enrichment analysis to be performed, run PathNet with Enrichment_Analysis=TRUE')
  }
  if (!Enrichment_Analysis && !Contextual_Analysis) {
    stop('You must specify that one or both of Enrichment_Analysis and Contextual_Analysis are TRUE')
  }
  
  # Construct DirectEvidence, get unique gene IDs, and calculate probabilities
  tmp.result <- construct.DirectEvidence(DirectEvidence_info, Column_DirectEvidence)
  DirectEvidence  <- tmp.result$DirectEvidence
  gene_ID         <- tmp.result$gene_ID
  
  # Adjust adjacency matrix
  Adjacency <- Adjacency[rownames(Adjacency) %in% gene_ID,rownames(Adjacency) %in% gene_ID]  
  observe <- rep(FALSE, nrow(Adjacency))
  names(observe) <- rownames(Adjacency)
  diag(Adjacency) <- 0
  
  # Prepare pathway
  pathway = pathway
  pathwayNames <- unique(pathway[,3])
  pathwayNames <- sort(pathwayNames)
  path_genes_all <- rep(0,1)
  
  # Processing variables
  ea.results  <- list(result=NULL, combined.evidence=NULL, rankwise=NULL)
  ctx.results <- list(conn_p_value=NULL, overlap=NULL)
  rankwise    <- NULL
  
  # **************************** Permutation based on all the genes ********************************************
  
  if (Enrichment_Analysis == TRUE)
  {
    ea.results <- do.enrichment.analysis(DirectEvidence, Adjacency, observe, gene_ID,
                                         pathway, pathwayNames, n_perm, threshold)
    # Update rankwise result potentially used during contextual analysis
    rankwise <- ea.results$rankwise
  }
  #*******************************************Calculation of pathway connectivity*******************************
  if (Contextual_Analysis == TRUE)
  {
    message('Starting contextual analysis...')
    ctx.results <- do.contextual.analysis(DirectEvidence, Adjacency, gene_ID, pathway, 
                                          pathwayNames, n_perm, use_sig_pathways, ea.results)
  }
  
  # Return objects
  list(enrichment_results=ea.results$result,
       enrichment_combined_evidence=ea.results$combined.evidence,
       conn_p_value=ctx.results$conn_p_value,
       pathway_overlap=ctx.results$overlap)
}


# Constructs initial Direct Evidence structures and determines which gene IDs to use
construct.DirectEvidence <- function(DirectEvidence_info, Column_DirectEvidence) {
  gene_ID <- as.numeric(DirectEvidence_info[, 1])
  message("Number of genes present in the DirectEvidence File:")
  message(length(gene_ID))
  
  DirectEvidence <- as.numeric(DirectEvidence_info[, Column_DirectEvidence])
  unique_gene_ids <- unique(gene_ID)
  
  # If multiple probes are present the probe with min p-value is selected
  if (length(gene_ID) != length(unique_gene_ids))
  {
    message("Finding unique genes...")
    temp_DirectEvidence <- rep(0, length(unique_gene_ids))    
    for (i in 1:length(unique_gene_ids))
    {
      temp_DirectEvidence[i] <- max(DirectEvidence[gene_ID %in% unique_gene_ids[i]])
    }
    gene_ID <- unique_gene_ids
    DirectEvidence <- as.matrix(temp_DirectEvidence)
    rownames(DirectEvidence) <- as.matrix(gene_ID)
    
    message(paste("Unique number of genes in DirectEvidence File:",length(gene_ID)))
    
    # TODO: Remove
    write.table(DirectEvidence, "Unique_gene_log10_pVal.txt", row.names = TRUE, sep = "\t")
    rm(unique_gene_ids,temp_DirectEvidence)
  }
  list(DirectEvidence=DirectEvidence, gene_ID=gene_ID)
}
