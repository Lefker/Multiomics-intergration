# Quantile_normalization
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

# Normalization function between and within columns
norm_fun <- function(df){
  dfl <- log2(df+1)
  df1 <- scale(dfl, center = TRUE, scale = TRUE)
  df2 <- scale(t(df1), center = TRUE, scale = TRUE)
  df3 <- t(df2)
  df4 <- quantile_normalisation(df3)
  return(df4)
}

# Fold Change
FoldChange_fun <- function(df1, df2, rnames){
  FCn = apply(df1,1,mean)
  FCc = apply(df2,1,mean)
  FC <- FCn - FCc
  FC <- as.matrix(FC)
  rownames(FC) <- rownames(rnames)
  colnames(FC) <- c("FoldChange")
  return(FC)
}

