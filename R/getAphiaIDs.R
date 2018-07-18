#' @title getAphiaIDs
#' @description This function contains all of the logic for passing the data off
#' to the various webservices that provide APHIAID
#' @param mysteryAPHIAID - df of species for which the APHIAID is either unknown or not
#' definitive
#' @param doSci - flag indicating whether scientific names should be used when 
#' checking the services
#' @param doComm - flag indicating whether common names should be used when 
#' checking the services
#' @param masterList - the (semi) original list of species provided.
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
getAphiaIDs<-function(mysteryAPHIAID = NULL, doSci=T, doComm=T, masterList = NULL){
  if (doSci & (nrow(mysteryAPHIAID)>0))  {
    recs_sci = unique(mysteryAPHIAID[!is.na(mysteryAPHIAID$SCI_COL_CLN),"SCI_COL_CLN"])
    if (length(recs_sci)>0) {
      cat(paste0("\ttaxize > scientific names\n"), file = "getTaxaIDs.log", append = TRUE)
      sci =   chk_taxize(recs_sci, searchtype = 'scientific')   
      sci = merge(mysteryAPHIAID[,-which(colnames(mysteryAPHIAID) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], sci, all.x=TRUE, by="SCI_COL_CLN")
      sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]," (SCIENTIFIC)")
      if (nrow(sci)>0) {
        defCheck = assignDefinitive(df = sci, masterList = masterList)
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        if (!exists("definitiveAPHIAID")){
          definitiveAPHIAID<-newdefinitive
        }else{
          definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
        rm(defCheck)
        rm(newdefinitive)
      }
    }

  }
  if (doComm & (nrow(mysteryAPHIAID)>0)) {
    recs_comm= unique(mysteryAPHIAID[!is.na(mysteryAPHIAID$COMM_COL_CLN),"COMM_COL_CLN"])
    if (length(recs_comm)>0) {
      cat(paste0("\ttaxize > common names\n"), file = "getTaxaIDs.log", append = TRUE)
      comm =   chk_taxize(recs_comm, searchtype = 'common')
    }
    
    comm = merge(mysteryAPHIAID[,-which(colnames(mysteryAPHIAID) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], comm, all.x=TRUE, by="COMM_COL_CLN")
    comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]," (COMMON)")
    if (nrow(comm)>0) {
      defCheck = assignDefinitive(df = comm, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryAPHIAID= defCheck[[2]]
      if (!exists("definitiveAPHIAID")){
        definitiveAPHIAID<-newdefinitive
      }else{
        definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
    }
  }
  if (doSci  & (nrow(mysteryAPHIAID)>0))  {
    recs_sci = unique(mysteryAPHIAID[!is.na(mysteryAPHIAID$SCI_COL_CLN),"SCI_COL_CLN"])
    if (length(recs_sci)>0) {
      cat(paste0("\tworrms > scientific names\n"), file = "getTaxaIDs.log", append = TRUE)
      sci =   chk_worrms(recs_sci, searchtype = 'scientific')
    }
    sci = merge(mysteryAPHIAID[,-which(colnames(mysteryAPHIAID) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], sci, all.x=TRUE, by="SCI_COL_CLN")
    sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]," (SCIENTIFIC)")
    if (nrow(sci)>0) {
      defCheck = assignDefinitive(df = sci, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryAPHIAID= defCheck[[2]]
      if (!exists("definitiveAPHIAID")){
        definitiveAPHIAID<-newdefinitive
      }else{
        definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      rm(defCheck)
      rm(newdefinitive)
    }
  }
  if (doComm & (nrow(mysteryAPHIAID)>0)) {    
    recs_comm= unique(mysteryAPHIAID[!is.na(mysteryAPHIAID$COMM_COL_CLN),"COMM_COL_CLN"])
    if (length(recs_comm)>0) {
      cat(paste0("\tworrms > common names\n"), file = "getTaxaIDs.log", append = TRUE)
      comm =   chk_worrms(recs_comm, searchtype = 'common')
    }
    
    comm = merge(mysteryAPHIAID[,-which(colnames(mysteryAPHIAID) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], comm, all.x=TRUE, by="COMM_COL_CLN")
    comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]," (COMMON)")
    if (nrow(comm)>0) {
      defCheck = assignDefinitive(df = comm, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryAPHIAID= defCheck[[2]]
      if (!exists("definitiveAPHIAID")){
        definitiveAPHIAID<-newdefinitive
      }else{
        definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
    }
  }
  res = list(definitiveAPHIAID, mysteryAPHIAID)
  return(res)
}