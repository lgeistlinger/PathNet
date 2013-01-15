library(PathNetData)
data(brain_regions)
data(A)
data(pathway)

# tests row-wise refactoring works as expected for calculating Direct Evidence
test_enrich.DE_rowise <- function() {
  evidence_column <- 7
  br <- brain_regions[100:500,]   # subset for speed
  tmp <- PathNet:::construct.DirectEvidence(br, evidence_column)
  DirectEvidence  <- tmp$DirectEvidence
  checkEquals(401, length(DirectEvidence))
  gene_ID <- tmp$gene_ID
  
  # setup
  Adjacency <- A[rownames(A) %in% gene_ID,rownames(A) %in% gene_ID]
  checkEquals(c(113, 113), dim(Adjacency))
  observe <- rep(FALSE, nrow(Adjacency))
  names(observe) <- rownames(Adjacency)
  diag(Adjacency) <- 0
  
  # original loop
  DER <- rep(NA, nrow(Adjacency))
  for (i in 1:nrow(Adjacency))
  {
    DER[i] <- DirectEvidence[gene_ID %in% names(observe[i])]
  }
  # new sapply 
  DER2 <- sapply(1:nrow(Adjacency), function(i) { DirectEvidence[gene_ID %in% names(observe[i]) ]})
  checkEquals(DER, DER2)
  checkEquals(0.683263367, DER[1], tolerance=0.00001)
  checkEquals(0.750759287, DER[113], tolerance=0.00001)
}

# Test calculating the pathway links
test_pathway_links <- function() {
  pathwayNames <- sort(unique(pathway[,3]))
  links <- PathNet:::construct.pathway_links(pathway, pathwayNames)
  # Check length and spot check counts
  checkEquals(130, length(links))
  checkEquals(0, links[1])
  checkEquals(1162, links[80])
  checkEquals(851, links[130])
}

# Test main enrichment function
test_enrichment <- function() {
  evidence <- brain_regions[1:100,]  # subset for speed requirements
  #evidence <- brain_regions
  
  set.seed(123)
  z <- PathNet(Enrichment_Analysis = TRUE, 
          DirectEvidence_info = evidence, 
          Adjacency = A, 
          pathway = pathway, 
          Column_DirectEvidence = 7,
          n_perm = 1, threshold = 0.05)
  
  # Enrichment results
  results <- z$enrichment_results
  checkEquals(130, nrow(results))
  idx <- 10
  checkEquals('Insulin signaling pathway', as.character(results$Name[idx]))
  checkEquals(2, results$No_of_Genes[idx])
  checkEquals(1, results$Sig_Direct[idx])
  checkEquals(1, results$Sig_Combi[idx])
  checkEquals(0.4857143, results$p_Hyper[idx], tolerance=0.00001)
  checkEquals(1, results$p_Hyper_FWER[idx])
  checkEquals(0.4563265, results$p_PathNet[idx], tolerance=0.00001)
  checkEquals(1, results$p_PathNet_FWER[idx])
  
  idx <- 36
  checkEquals('Vascular smooth muscle contraction', as.character(results$Name[idx]))
  checkEquals(13, results$No_of_Genes[idx])
  checkEquals(3, results$Sig_Direct[idx])
  checkEquals(3, results$Sig_Combi[idx])
  checkEquals(0.7900371, results$p_Hyper[idx], tolerance=0.00001)
  checkEquals(1, results$p_Hyper_FWER[idx])
  checkEquals(0.7341656, results$p_PathNet[idx], tolerance=0.00001)
  checkEquals(1, results$p_PathNet_FWER[idx])
  
  # Combined and Indirect Evidence
  er <- z$enrichment_combined_evidence
  idx <- 23
  checkEquals(31, er[idx,1])
  checkEquals(0.0881113798, er[idx,2], tolerance=0.00001)
  checkTrue(is.na(er[idx,3]))
  checkEquals(0.088111380, er[idx,4], tolerance=0.00001)
}
