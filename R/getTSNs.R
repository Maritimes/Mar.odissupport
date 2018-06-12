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
  print("Finding TSNs")
  if (!is.null(knownAphias))  {
    print(paste0("::worrms::"))
    print(paste0("---Existing AphiaIDs---"))
    
    aphiaTSNRes =   chk_worrmsTSN(knownAphias)
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
    print(paste0("::ritis::"))
    print(paste0("---scientific names---"))
    
    cln = skipUselessRecs(mysteryTSN, "SCI_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      sci = cln[[1]]
    }else{
      sci =   chk_ritis(cln[[1]], "SCI_COL_CLN",searchtype = 'scientific')
    }
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
      if (nrow(cln[[2]])>0){
        mysteryTSN = rbind(mysteryTSN,cln[[2]])
        rm(cln)
      }
    }
  }
  if (doComm  & (nrow(mysteryTSN)>0))  {
    print(paste0("---common names---"))
    cln = skipUselessRecs(mysteryTSN, "COMM_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      comm = cln[[1]]
    }else{
      comm =   chk_ritis(cln[[1]], "COMM_COL_CLN",searchtype = 'common')
    }

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
      if (nrow(cln[[2]])>0){
        mysteryTSN = rbind(mysteryTSN,cln[[2]])
        rm(cln)
      }
    }
  }
  res = list(definitiveTSN, mysteryTSN)
  return(res)
}