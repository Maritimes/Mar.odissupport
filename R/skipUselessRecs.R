#' @title skipUselessRecs
#' @description This function is used to prevent records from being sent to
#' webservices if they do not contain anything that could be used to determine
#' a definitive code.
#' @param df This is a dataframe of species for which we have no definitive code
#' @param field - field in the df containing the value to check against the 
#' service
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
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