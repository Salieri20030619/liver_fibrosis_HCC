u_stat <- function(ranks_matrix, gene_idx, maxRank=1500, sparse=FALSE){
  
    len_sig <- length(gene_idx)
    ncells <- ncol(ranks_matrix)-1 
    
    if (len_sig <= 0) {
      return(rep(0, ncells))
    }
    present_idx <- gene_idx[gene_idx>0]
    missing_idx <- gene_idx[gene_idx<0]
  
    #Initialize rank_sum with missing genes
    rank_sum <- rep(length(missing_idx) * maxRank, ncells)
  
    #Subset ranks_matrix by rows (genes) using numeric indices
    if (length(present_idx) > 0) {
      rank_sub <- ranks_matrix[present_idx, -1, drop = FALSE]
      if (sparse==TRUE){
        rank_sub[rank_sub==0] <- maxRank
      }
      insig <- rank_sub >= maxRank
      rank_sub[insig] <- maxRank
    
      rank_sum = rank_sum + apply(rank_sub, 2, sum)
    }
  
    rank_sum_min <- len_sig*(len_sig + 1)/2
    ucell_score <- 1 - (rank_sum-rank_sum_min)/(len_sig*maxRank - rank_sum_min)
    return(ucell_score)
}
x<-read.table("/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/analysis_stage1_20261017_075157UTC/F4_UCell_Hallmark_v1/test_ranks.tsv",header=FALSE); x<-cbind(rn=seq_len(nrow(x)),x); write.table(u_stat(x,c(1,3,-1),maxRank=4),file="/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/analysis_stage1_20261017_075157UTC/F4_UCell_Hallmark_v1/R_expected.tsv",row.names=FALSE,col.names=FALSE)
