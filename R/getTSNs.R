#' @title getTSNs
#' @description This function contains all of the logic for passing the data off
#' to the various webservices that provide TSNs.
#' @param mysteryTSN - df of species for which the TSN is either unknown or not
#' definitive
#' @param doSci - flag indicating whether scientific names should be used when 
#' checking the services
#' @param doComm - flag indicating whether common names should be used when 
#' checking the services
#' @param knownAphias - one worrms check finds TSNs from AphiaIDs 
#' @param masterList - the (semi) original list of species provided.
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
getTSNs<-function(mysteryTSN = NULL, doSci=T, doComm=T, knownAphias = NULL, masterList = NULL){
  if (!is.null(knownAphias))  {
    recs_aphiaid = unique(knownAphias[!is.na(knownAphias$CODE),"CODE"])
    
    aphiaTSNRes =   chk_worrmsTSN(recs_aphiaid,knownAphias)
    colnames(knownAphias)[colnames(knownAphias)=='CODE']<-"APHIAID"
    aphiaTSNRes = merge(knownAphias[,c("ID","APHIAID","SCI_COL_CLN","COMM_COL_CLN")],aphiaTSNRes, by = "APHIAID")
    aphiaTSNRes$APHIAID<-NULL
   
    defCheck = assignDefinitive(df = aphiaTSNRes, masterList = masterList)
    newdefinitive= defCheck[[1]]
    #mysteryTSN= defCheck[[2]] #these are just the ones we couldn't find via an aphiaid
    mysteryTSN= rbind(mysteryTSN[!(mysteryTSN$ID %in% knownAphias$ID),],defCheck[[2]])
    if (!exists("definitiveTSN")){
      definitiveTSN<-newdefinitive
    }else{
      definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
    }
    rm(defCheck)
    rm(newdefinitive)
  }
  if (doSci  & (nrow(mysteryTSN)>0))  {

    recs_sci = unique(mysteryTSN[!is.na(mysteryTSN$SCI_COL_CLN),"SCI_COL_CLN"])
    if (length(recs_sci)>0) {
      sci =   chk_ritis(recs_sci, searchtype = 'scientific')
    }
    sci = merge(mysteryTSN[,-which(colnames(mysteryTSN) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], sci, all.x=TRUE, by.x="SCI_COL_CLN", by.y="joincol")
    sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(sci[!is.na(sci$SUGG_SPELLING),"SUGG_SPELLING"]," (SCIENTIFIC)")

    if (nrow(sci)>0) {
      defCheck = assignDefinitive(df = sci, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryTSN= defCheck[[2]]
      if (!exists("definitiveTSN")){
        definitiveTSN<-newdefinitive
      }else{
        definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      rm(defCheck)
      rm(newdefinitive)
    }
  }
  if (doComm  & (nrow(mysteryTSN)>0))  {
    recs_comm = unique(mysteryTSN[!is.na(mysteryTSN$COMM_COL_CLN),"COMM_COL_CLN"])
    if (length(recs_comm)>0) {
      comm =   chk_ritis(recs_comm, searchtype = 'common')
    }
    comm = merge(mysteryTSN[,-which(colnames(mysteryTSN) %in% c("CODE","CODE_SVC","CODE_TYPE","CODE_SRC","CODE_DEFINITIVE","SUGG_SPELLING"))], comm, all.x=TRUE, by.x="COMM_COL_CLN", by.y="joincol")
    comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]<-paste0(comm[!is.na(comm$SUGG_SPELLING),"SUGG_SPELLING"]," (COMMON)")

    if (nrow(comm)>0) {
      defCheck = assignDefinitive(df = comm, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryTSN= defCheck[[2]]
      if (!exists("definitiveTSN")){
        definitiveTSN<-newdefinitive
      }else{
        definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      rm(defCheck)
      rm(newdefinitive)
    }
  }
  res = list(definitiveTSN, mysteryTSN)
  return(res)
}