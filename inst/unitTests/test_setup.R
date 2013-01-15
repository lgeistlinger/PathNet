library(PathNetData)

test_construct.DirectEvidence <- function() {
  # Tests constructing direct evidence and gene IDs
  # Uses subset of brain_region data for speed
  data(brain_regions)
  br <- brain_regions
  
  tmp <- PathNet:::construct.DirectEvidence(br, 7)
  #tmp <- construct.DirectEvidence(br, 7)
  de  <- tmp$DirectEvidence
  gid <- tmp$gene_ID
  
  checkEquals(1.44304379, de[1])
  checkEquals(0.68326337, de[100])
  checkEquals(20077, length(gid))
  checkEquals(140:143, gid[96:99])
  checkEquals(146, gid[100])
}

# Test that we must specify some work to perform
test_enrichment_or_contextual <- function() {
  evidence <- brain_regions
  set.seed(123)
  tryCatch(PathNet(Enrichment_Analysis = FALSE,
                   Contextual_Analysis = FALSE,
                   DirectEvidence_info = evidence, 
                   Adjacency = A, 
                   pathway = pathway, 
                   Column_DirectEvidence = 7,
                   n_perm = 1, threshold = 0.05,
                   use_sig_pathways=TRUE),
           error = function(e) print(e$message))
}

# Test that we handle invalid run settings correctly
test_contextual_sig_pathways_requires_enrichment <- function() {
  evidence <- brain_regions
  set.seed(123)
  tryCatch(PathNet(Enrichment_Analysis = FALSE,
                   Contextual_Analysis = TRUE,
                   DirectEvidence_info = evidence, 
                   Adjacency = A, 
                   pathway = pathway, 
                   Column_DirectEvidence = 7,
                   n_perm = 1, threshold = 0.05,
                   use_sig_pathways=TRUE),
           error = function(e) print(e$message))
}
