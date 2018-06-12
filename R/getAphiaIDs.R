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
  print("Finding AphiaIDs")
  print(paste0("::Taxize::"))
  if (doSci & (nrow(mysteryAPHIAID)>0))  {
    print(paste0("---scientific names---"))
    
    cln = skipUselessRecs(mysteryAPHIAID, "SCI_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      sci = cln[[1]]
    }else{
      sci =   chk_taxize(cln[[1]], "SCI_COL_CLN",searchtype = 'scientific')
    }
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
      if (nrow(cln[[2]])>0){
        mysteryAPHIAID = rbind(mysteryAPHIAID,cln[[2]])
        rm(cln)
      }
    }
  }
  
  if (doComm & (nrow(mysteryAPHIAID)>0)) {
    print(paste0("---common names---"))
    
    cln = skipUselessRecs(mysteryAPHIAID, "COMM_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      comm = cln[[1]]
    }else{
      comm =   chk_taxize(cln[[1]], "COMM_COL_CLN",searchtype = 'common')
    }

    if (nrow(comm)>0) {
      defCheck = assignDefinitive(df = comm, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryAPHIAID= defCheck[[2]]
      if (!exists("definitiveAPHIAID")){
        definitiveAPHIAID<-newdefinitive
      }else{
        definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      rm(defCheck)
      rm(newdefinitive)
      if (nrow(cln[[2]])>0){
        mysteryAPHIAID = rbind(mysteryAPHIAID,cln[[2]])
        rm(cln)
      }
    }
  }
  
  print(paste0("::worrms::"))
  #worrms can generate multiple results/searchterm
  if (doSci  & (nrow(mysteryAPHIAID)>0))  {
    print(paste0("---scientific names---"))
    
    cln = skipUselessRecs(mysteryAPHIAID, "SCI_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      sci = cln[[1]]
    }else{
      sci =   chk_worrms(cln[[1]], "SCI_COL_CLN", searchtype = 'scientific')
    }
    
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
      if (nrow(cln[[2]])>0){
        mysteryAPHIAID = rbind(mysteryAPHIAID,cln[[2]])
        rm(cln)
      }
    }
  }
  if (doComm & (nrow(mysteryAPHIAID)>0)) {
    print(paste0("---common names---"))
    cln = skipUselessRecs(mysteryAPHIAID, "COMM_COL_CLN")
    if (nrow(cln[[1]])<1){
      print("No valid values to check - skipping check")
      comm = cln[[1]]
    }else{
      comm =   chk_worrms(cln[[1]], "COMM_COL_CLN", searchtype = 'common')
    }
    if (nrow(comm)>0) {
      defCheck = assignDefinitive(df = comm, masterList = masterList)
      newdefinitive= defCheck[[1]]
      mysteryAPHIAID= defCheck[[2]]
      if (!exists("definitiveAPHIAID")){
        definitiveAPHIAID<-newdefinitive
      }else{
        definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      rm(defCheck)
      rm(newdefinitive)
      if (nrow(cln[[2]])>0){
        mysteryAPHIAID = rbind(mysteryAPHIAID,cln[[2]])
        rm(cln)
      }
    }
  }
  res = list(definitiveAPHIAID, mysteryAPHIAID)
  return(res)
}