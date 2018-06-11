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
    sci =   chk_taxize(mysteryAPHIAID, "SCI_COL_CLN",searchtype = 'scientific')
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
    print(paste0("---common names---"))
    comm = chk_taxize(mysteryAPHIAID, "COMM_COL_CLN", searchtype = 'common')
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
    }
  }
  
  print(paste0("::worrms::"))
  #worrms can generate multiple results/searchterm
  if (doSci  & (nrow(mysteryAPHIAID)>0))  {
    print(paste0("---scientific names---"))
    sci =   chk_worrms(mysteryAPHIAID, "SCI_COL_CLN", searchtype = 'scientific')
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
    print(paste0("---common names---"))
    comm =   chk_worrms(mysteryAPHIAID, "COMM_COL_CLN", searchtype = 'common')
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
    }
  }
  res = list(definitiveAPHIAID, mysteryAPHIAID)
  return(res)
}