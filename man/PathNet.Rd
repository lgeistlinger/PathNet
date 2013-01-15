\name{PathNet}
\alias{PathNet}
\title{Enrichment and contextual analysis of KEGG pathways}
\description{
Pathway analysis using Network information (PathNet) algorithm, described in Dutta et al., is implemented in this package. PathNet uses topological information present in pathways and differential expression levels of genes (obtained from microarray experiment) to identify pathways that are 1) significantly enriched and 2) contextually associated with each other. In enrichment analysis, PathNet considers both the differential expression of genes and their pathway neighbors to strengthen the evidence that a pathway is implicated in the biological conditions tested for in the experiment. In addition, PathNet also identifies contextual association between pathway pairs based on connectivity of differentially expressed genes between two pathways.}

\usage{PathNet(Enrichment_Analysis = TRUE, Contextual_Analysis = FALSE,
               DirectEvidence_info, Column_DirectEvidence = 7, 
               Adjacency, pathway = pathway, n_perm = 2000,
               threshold = 0.05, use_sig_pathways  = FALSE)}
\arguments{
  \item{Enrichment_Analysis}{If set to TRUE, enrichment analysis will be carried out}
  \item{Contextual_Analysis}{If set to TRUE, contextual analysis will be carried out}
  \item{DirectEvidence_info}{Data containing direct evidence information. The negative log10 transformed p-value of the significance of association is used as direct evidence.}
  \item{Column_DirectEvidence}{Column of the data matrix to be used}
  \item{Adjacency}{Adjacency matrix}
  \item{pathway}{Pathway data}
  \item{n_perm}{Number of permutations}
  \item{threshold}{Threshold of significance}
  \item{use_sig_pathways}{A logical variable whether only significance pathways will be used for contextual analysis.  If this is TRUE then the Enrichment_Analysis parameter must also be set to TRUE.}
}

\value{
  The PathNet functions returns a list of values reflecting the results of performing pathway
  enrichment and/or contextual analysis.  The items returned in the list are:
  
  \item{enrichment_results}{Significance levels of enrichment for each pathway from PathNet and the hypergeometric test. These levels are in a dataframe which contains gene ID, direct, indirect, and combined evidence levels of genes present in the microarray data. Indirect evidences are calculated only for the genes that are present in the KEGG pathway and have at least one edge. For the rest of the genes, indirect evidences are ``NA'' and the combined evidences are replaced by the direct evidences.}
  
  \item{enrichment_combined_evidence}{A data frame containing gene ID, direct, indirect, and combined evidence levels of genes present in the microarray data.
  The results are sorted in increasing order based on the p_PathNet_FWER value. Results include pathway name (Name), number of genes of the pathway present in the microarray data (No_of_Genes), number of genes significant from direct evidence (Sig_Direct), number of genes significant from combined evidence (Sig_Combi), significance of enrichment from hypergeometric test (p_Hyper), family wise error rate correction of p_Hyper (p_Hyper_FWER), significance of enrichment from PathNet (p_PathNet), and family wise error rate correction of p_PathNet (p_PathNet_FWER).}
  
  \item{conn_p_value}{Contains the results of the contextual association in the form of a square matrix, where the number of rows (and columns) is equal to the number of pathways. The element at the ith row and jth column denotes the significance of contextual association of pathway i with pathway j.}
  
  \item{pathway_overlap}{The statistical significance of the overlapping genes between all pathway pairs in the form of a square matrix. This information is only based on the KEGG database and is not dependent on gene expression data. The hypergeometric test is used to estimate if the observed overlap is statistically significant. If use_sig_pathways parameter is set to TRUE, both contextual association and overlap analysis results are sorted in decreasing order of pathway enrichment. Otherwise, the pathways appear in alphabetical order of pathway names.}
}

\details{See cited document for more details}

\references{
PathNet: A tool for pathway analysis using topological information.
Dutta B, Wallqvist A, and Reifman J.
Source Code for Biology and Medicine 2012 Sep 24;7(1):10.
}
\keyword{Pathway}

\examples{
library(PathNetData)
data(brain_regions)
data(pathway)
data(A)

# Perform an example enrichment analysis using  a
# subset of data and a small number of permutations
# for demonstration purposes.
results <- PathNet(Enrichment_Analysis = TRUE, Contextual_Analysis = FALSE, 
        DirectEvidence_info = brain_regions[1:100,], Column_DirectEvidence = 7,
        Adjacency = A , pathway = pathway, n_perm = 10, threshold = 0.05,
        use_sig_pathways  = FALSE)
        
# display part of results
head(results$enrichment_results)
head(results$enrichment_combined_evidence)
}
