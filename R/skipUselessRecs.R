skipUselessRecs <-function(df = NULL, field = NULL){
  skips = df[0,]
  if (nrow(df[is.na(df[,field]),])>0)skips = rbind(skips, df[is.na(df[,field]),])
  if (nrow(df[nchar(df[,field])<4,])>0)skips = rbind(skips, df[nchar(df[,field])<4,])
  mystery = df[!(df$ID %in% skips$ID),]

  if (nrow(skips)>0) {
    skips = skips[!is.na(skips$ID),]
  }
  res=list(mystery, skips)
  return(res)
}