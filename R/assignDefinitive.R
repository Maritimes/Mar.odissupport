#' @title assignDefinitive
#' @description This function takes a dataframe of results from a webservice, 
#' and determines whether the newly found values are definitive, or if they 
#' should be sent to additional web services in an attempt to get a more 
#' conclusive result.
#' If multiple services return the same code, the results are aggregated so that 
#' all of the contributing services will be captured (e.g. this code was the 
#' result of taxize AND worrms)
#' @param df This is a dataframe of species for which we have no definitive code
#' @param masterList  This is the originally submitted list of species (with an 
#' internally assigned ID)
#' @importFrom stats aggregate
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
assignDefinitive <- function(df = NULL, masterList = NULL){
  fields = names(masterList)
  unknown = df[!(df$CODE_DEFINITIVE %in% TRUE), ]
  known = merge(masterList[masterList$ID %in% df[df$CODE_DEFINITIVE %in% TRUE,"ID"],fields],
                df[df$CODE_DEFINITIVE %in% TRUE,c("ID","CODE","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")], by.x = "ID", by.y = "ID", all.x=T)
  
  unknown = merge(masterList, unknown[, c("ID","CODE","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")], by.x = "ID", by.y = "ID", all.y=T) 
  unknown = unique(unknown)
  #remove known values from mystery
  unknown = unknown[!(paste0(unknown$ID,unknown$CODE_TYPE) %in% paste0(known$ID, known$CODE_TYPE)),]
  #remove duplicates
  probs = unique(unknown[,c("CODE","ID","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")])
  if (nrow(probs) < nrow(unknown)){
    browser()
    unknownAgg = aggregate(by=unknown[c("CODE","ID","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")], x = unknown[c("CODE_SVC")], paste, collapse = ",")
    unknownAgg = merge(masterList,unknownAgg, all.y = T)
    unknownFinal = unique(rbind(unknownAgg,unknown[!(unknown$ID %in% unknownAgg$ID),]))
  }else{
    unknownFinal = unique(unknown)
  }
  resAss = list(known, unknownFinal)
  return(resAss)
}