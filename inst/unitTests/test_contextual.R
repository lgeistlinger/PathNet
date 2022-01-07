library(PathNetData)
data(brain_regions)
data(A)
data(pathway)


# Test determining common path
test_determine.common.path <- function() {
  # Setup variables
  pathwayNames <- unique(pathway[,3])
  pathwayNames_significant <- sort(pathwayNames)
  tmp <- PathNet:::construct.DirectEvidence(brain_regions, 7)
  DirectEvidence  <- tmp$DirectEvidence
  gene_ID <- tmp$gene_ID
  Adjacency <- A
  Adjacency <- Adjacency[rownames(Adjacency) %in% gene_ID,rownames(Adjacency) %in% gene_ID]  
  
  # Test function
  common_path_results <- PathNet:::determine.common.path(pathwayNames_significant[1],
                                               pathway, DirectEvidence, 
                                               Adjacency, gene_ID)
  path_genes1                 <- common_path_results$path_genes1
  path_genes_common1          <- common_path_results$path_genes_common1
  DirectEvidence_path_common1 <- common_path_results$DirectEvidence_path_common1
  
  checkEquals(c(8647, 64241, 64240, 9429), path_genes1[1:4])
  checkEquals(c(1080, 6833, 6891, 6890), as.integer(path_genes_common1))
}

# Tests main contextual analysis operation WITHOUT enrichment analysis
test_contextual_no_enrichment <- function() {
  evidence <- brain_regions[1:100,]   # subset for speed
  set.seed(123)
  z <- PathNet(Enrichment_Analysis = FALSE,
               Contextual_Analysis = TRUE,
               DirectEvidence_info = evidence, 
               Adjacency = A, 
               pathway = pathway, 
               Column_DirectEvidence = 7,
               n_perm = 5, threshold = 0.05)
  checkEquals(NULL,z$enrichment_results)
  checkEquals(NULL,z$enrichment.combined.evidence)
  
  cp <- z$conn_p_value
  checkEquals(c(0,1,1), cp[1,1:3], checkNames=FALSE)
  
  po <- z$pathway_overlap
  idx <- 5
  checkEquals(1, po[idx,1])
  checkEquals(3.291294e-12, po[idx,2], tolerance=0.00001)
  checkEquals(0.0332455496, po[idx,3], tolerance=0.00001)
  checkEquals(2.484184e-02, po[idx,4], tolerance=0.00001)
  checkEquals(0, po[idx,5], tolerance=0.00001)
}

# Test construct sig pathways function
test_construct.significant.pathways_nosig <- function() {
  # use_sig_pathways, pathwayNames, ea.results
  pathwayNames <- unique(pathway[,3])
  pathwayNames <- sort(pathwayNames)
  # first argument FALSE means we just sort
  sigpaths <- PathNet:::construct.significant.pathways(FALSE, pathwayNames, NULL)
  checkEquals(130, length(pathwayNames))
  checkEquals(130, length(sigpaths))
  checkEquals('ABC transporters', as.character(sigpaths[1]))
  checkEquals('Wnt signaling pathway', as.character(sigpaths[130]))
  checkEquals(as.character(pathwayNames), as.character(sigpaths))
}

# Full test contextual analysis using signficant pathways
#test_contextual_with_sig_pathways <- function() {
#  evidence <- brain_regions[1:5000,]   # limit number of regions for speed
#  set.seed(123)
#
#  z <- PathNet(Enrichment_Analysis = TRUE,
#               Contextual_Analysis = TRUE,
#               DirectEvidence_info = evidence, 
#               Adjacency = A, 
#               pathway = pathway, 
#               Column_DirectEvidence = 7,
#               n_perm = 10, threshold = 0.05,
#               use_sig_pathways=TRUE)
#  
#  # Check p-values
#  cp <- z$conn_p_value
#  checkEquals(c(0.0, 0.0, 0.0, 0.2, 1.0), cp[2,1:5], checkNames=FALSE)
#  
#  # Pathway overlap
#  po <- z$pathway_overlap
#  idx <- 2
#  checkEquals('Acute myeloid leukemia', rownames(po)[idx])
#  checkEquals('Adherens junction', colnames(po)[1])
#  checkEquals( 0.0004628949, po[idx,1], tolerance=0.00001)
#  checkEquals( 0, po[idx,2])
#  checkEquals( 1.284468e-08, po[idx,3], tolerance=0.00001)
#  checkEquals( 3.291294e-12, po[idx,4], tolerance=0.00001)
#  checkEquals( 1, po[idx,5])
#}
