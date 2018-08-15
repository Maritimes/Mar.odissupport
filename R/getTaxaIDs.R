#' @title getTaxaIDs
#' @description This function attempts to determine a definitive AphiaID and TSN
#' from various web services via some existing R packages including worrms, 
#' taxize and ritis.
#' 
#' It takes a species list with Scientific and Common names and performs checks
#' in the following order:
#' \enumerate{
#'   \item Taxize (for AphiaID): Scientific name
#'   \item Worrms (for AphiaID): Scientific name
#'   \item Taxize (for AphiaID): Common name
#'   \item Worrms (for AphiaID): Common name
#'   \item Worrms (for TSN): <previously found> AphiaIDs
#'   \item Ritis (for TSN): Scientific name
#'   \item Ritis (for TSN): Common name
#'   }
#'   
#' In addition to all of the original fields, the returned data
#' will include: 
#' \enumerate{
#'   \item APHIAID
#'   \item APHIAID_SRC - what data was used to find the Aphiaid value (e.g. 
#'   scientific name, common name)
#'   \item APHIAID_SVC - which service(s) provided the Aphiaid value (e.g. 
#'   taxize, ritis, worrms)
#'   \item APHIAID_DEFINITIVE - TRUE indicates a single, confident match with a 
#'   service, FALSE indicates that either several potential matches were found, 
#'   or that the matches were recognized as authoritative by the service
#'   \item APHIAID_SPELLING - alternative spellings suggested for the 
#'   APHIAID_SRC suggested by the service that found the APHIAID
#'   \item TSN 
#'   \item TSN_SRC - what data was used to find the TSN value (e.g. 
#'   scientific name, common name, APHIAID)
#'   \item TSN_SVC - which service(s) provided the TSN value (e.g. 
#'   taxize, ritis, worrms)
#'   \item TSN_definitive - TRUE indicates a single, confident match with a 
#'   service, FALSE indicates that either several potential matches were found, 
#'   or that the matches were recognized as authoritative by the service
#'   \item TSN_SPELLING - alternative spellings suggested for the 
#'   TSN_SRC suggested by the service that found the TSN
#'   }
#' @param spec_list the dataframe containing information to be decoded to TSN and
#' aphiaIDs
#' @param sci_col the name of the column of the dataframe containing the
#' scientific names
#' @param comm_col the name of the column of the dataframe containing the
#' common names
#' @param sci_Filts default is \code{NULL} -  a vector of regex values that you 
#' might want to filter out of your scientific names prior to sending them to a 
#' web service.  For example, some Maritimes names inlude "(NS)", which will 
#' prevent services from finding matches, By adding "\\(NS\\)" (escaping the 
#' brackets and periods with slashes), we ensure that the results will be as 
#' clean as possible prior to searching. 
#' @param comm_Filts default is \code{NULL} -  a vector of regex values that you 
#' might want to filter out of your common names prior to sending them to a 
#' web service.  
#' By default, the following values will be filtered from both \code{sci_col} 
#' and \code{comm_col}:
#' \itemize{
#' \item \\(.*?\\) -removes thing in brackets (e.g "(NS)")
#' \item \\b[a-zA-Z]{1,2}\\. - removes one or two letter blocks of text 
#' (potentially followed by a period) (e.g "SP.")
#' \item \\,\\s?(SMALL|LARGE) - removes instances like ",SMALL"
#' \item UNIDENTIFIED - removes the word "UNIDENTIFIED"
#' \item UNID\\. - removes the word "UNID."
#' \item EGGS - removes the word "EGSS"
#' \item PURSE\\s" - removes the word "PURSE"
#' }
#' @param codes These are the codes you would like to determine.  \code{TSN} and
#' \code{APHIAID} are the only valid entries.
#' @param debug default is \code{FALSE}  This just ensure that the log file is 
#' overwritten rather than making many new ones.
#' @family speciesCodes
#' @importFrom stats setNames
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
getTaxaIDs <- function(spec_list = NULL,
                       sci_col = NULL,
                       comm_col = NULL,
                       sci_Filts = NULL,
                       comm_Filts = NULL,
                       codes = c("APHIAID", "TSN"),
                       debug = F) {
  ts = ifelse(debug, "000", format(Sys.time(), "%Y%m%d_%H%M"))
  
  logName = paste0("getTaxaIDs_",ts,".log")
  log_con <- file(logName)
  start = Sys.time()
  cat(paste0("Initiated at ",format(start, "%Y-%m-%d %H:%M"),"\n"), file = log_con) 
  cat(paste0("A log file has been generated at ",file.path(getwd(),logName),"\nIn the event of failure, this file will reveal the last successful operation\n"))
  fields = c("ID")
  
  
  
  if (!is.null(sci_col)){
    spec_list$SCI_COL_CLN = NA
    fields = c(fields, "SCI_COL_CLN")
    doSci = T  
  } else{
    doSci = F
  }
  if (!is.null(comm_col)){
    spec_list$COMM_COL_CLN = NA
    fields = c(fields, "COMM_COL_CLN")
    doComm = T  
  } else{
    doComm = F
  }
  
  spec_list = cleanPrepareSpecList(spec_list, sci_col, comm_col, sci_Filts, comm_Filts)
  cat(paste0("Input data successfully filtered\n"), file = logName, append = TRUE)
  
  definitive = spec_list[, fields]
  
  #add some columns to be populated later
  cols = c("CODE","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")
  definitive = cbind(definitive, setNames( lapply(cols, function(x) x=NA), cols) )
  definitive$CODE_DEFINITIVE = NA
  
  if ("SCI_COL_CLN" %in% names(spec_list) & "COMM_COL_CLN" %in% names(spec_list)){
    #ensure that each record has at least one valid value
    definitive = definitive[!is.na(definitive$SCI_COL_CLN) | !is.na(definitive$COMM_COL_CLN) ,]
  } else if ("SCI_COL_CLN" %in% names(spec_list)){
    definitive = definitive[!is.na(definitive$SCI_COL_CLN) ,]
  } else if ("COMM_COL_CLN" %in% names(spec_list)){
    definitive = definitive[!is.na(definitive$COMM_COL_CLN) ,]
  } 
  mysteryAPHIAID<-mysteryTSN<-definitive
  knownAPHIAID<-uncertainAPHIAID<-knownTSN<-uncertainTSN<-data.frame()
  rm(definitive)
  
  #approach -- Scientific names are the most likely source of a definitive match
  #so we look for them first (for both AphiaIDs and TSNs)
  # Check scientific names for AphiaIDs using taxize, then worrms
  # if TSNs have been requested, use the aphiaids to look for them
  
  TSNs_checked = c()
  APHIAIDs_checked = c()
  
  if(doSci & "APHIAID" %in% codes) {
    APHIAID_SCI =   getCodes(mysteryAPHIAID, definitiveAPHIAID = NULL, APHIAIDs_checked, mysteryTSN, definitiveTSN = NULL, TSNs_checked, "scientific", "APHIAID", codes, spec_list, logName)
    definitiveAPHIAID = APHIAID_SCI[[2]]
    definitiveTSN = APHIAID_SCI[[4]]
    mysteryAPHIAID = APHIAID_SCI[[1]]
    mysteryTSN = APHIAID_SCI[[3]]
    APHIAIDs_checked = APHIAID_SCI[[5]]
    TSNs_checked = APHIAID_SCI[[6]]
    mysteryAPHIAID = mysteryAPHIAID[!(mysteryAPHIAID$ID %in% definitiveAPHIAID$ID),]
    mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% definitiveTSN$ID),]
    # print(mysteryAPHIAID)
    # print(definitiveAPHIAID)
     #print(mysteryTSN)
     #print(definitiveTSN)
    #print("####")
  }
  if(doSci & "TSN" %in% codes){
    TSN_SCI =           getCodes(mysteryAPHIAID, definitiveAPHIAID, APHIAIDs_checked, mysteryTSN, definitiveTSN, TSNs_checked, "scientific", "TSN", codes, spec_list, logName)
    definitiveAPHIAID = TSN_SCI[[2]]
    definitiveTSN = TSN_SCI[[4]]
    mysteryAPHIAID = TSN_SCI[[1]]
    mysteryTSN = TSN_SCI[[3]]
    APHIAIDs_checked = TSN_SCI[[5]]
    TSNs_checked = TSN_SCI[[6]]
    mysteryAPHIAID = mysteryAPHIAID[!(mysteryAPHIAID$ID %in% definitiveAPHIAID$ID),]
    mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% definitiveTSN$ID),]
    # print(mysteryAPHIAID)
    # print(definitiveAPHIAID)
     #print(mysteryTSN)
     #print(definitiveTSN)
    #print("####")
  }
  if(doComm & "APHIAID" %in% codes) {
    APHIAID_COMM = getCodes(mysteryAPHIAID, definitiveAPHIAID, APHIAIDs_checked, mysteryTSN, definitiveTSN, TSNs_checked, "common", "APHIAID", codes, spec_list, logName)
    definitiveAPHIAID = APHIAID_COMM[[2]]
    definitiveTSN = APHIAID_COMM[[4]]
    mysteryAPHIAID = APHIAID_COMM[[1]]
    mysteryTSN = APHIAID_COMM[[3]]
    APHIAIDs_checked = APHIAID_COMM[[5]]
    TSNs_checked = APHIAID_COMM[[6]]
    mysteryAPHIAID = mysteryAPHIAID[!(mysteryAPHIAID$ID %in% definitiveAPHIAID$ID),]
    mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% definitiveTSN$ID),]
    # print(mysteryAPHIAID)
    # print(definitiveAPHIAID)
    #print(mysteryTSN)
     #print(definitiveTSN)
    #print("####")
  }
  if(doComm & "TSN" %in% codes) {
    TSN_COMM =         getCodes(mysteryAPHIAID, definitiveAPHIAID, APHIAIDs_checked, mysteryTSN, definitiveTSN, TSNs_checked, "common", "TSN", codes, spec_list, logName)
    definitiveAPHIAID = TSN_COMM[[2]]
    definitiveTSN = TSN_COMM[[4]]
    mysteryAPHIAID = TSN_COMM[[1]]
    mysteryTSN = TSN_COMM[[3]]
    APHIAIDs_checked = TSN_COMM[[5]]
    TSNs_checked = TSN_COMM[[6]]
    mysteryAPHIAID = mysteryAPHIAID[!(mysteryAPHIAID$ID %in% definitiveAPHIAID$ID),]
    mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% definitiveTSN$ID),]
    # print(mysteryAPHIAID)
    # print(definitiveAPHIAID)
     #print(mysteryTSN)
     #print(definitiveTSN)
    #print("####")
  }
  
  #This was going to be an attempt to use suggested spellings from services to  
  #re-initiate some of the checks.  I decided it's getting too automated and 
  #it will likely increase false positives.  For example, "STONE" will find
  #many valid values, like "STONEFLIES", "STONEFISH", etc, and it is difficult
  #to automate their correct assignment.  For now, such results remain as 
  #multi-results.
  
  # if (nrow(mysteryAPHIAID[!is.na(mysteryAPHIAID$SUGG_SPELLING),])>0 |
  #     nrow(mysteryTSN[!is.na(mysteryTSN$SUGG_SPELLING),])>0){
  #   start_spell = Sys.time()
  #   cat(paste0("Reviewing discovered alternative spellings to see if anything new can be found\n"), file = logName, append = TRUE)
  #   spell_APHIAID <- mysteryAPHIAID[!is.na(mysteryAPHIAID$SUGG_SPELLING),]
  #   if (nrow(spell_APHIAID)>0){
  #     #browser()
  #     # for (i in 1:nrow(spell_APHIAID)){
  #     #   
  #     # }
  #   }
  #   spell_TSN <- mysteryTSN[!is.na(mysteryTSN$SUGG_SPELLING),]
  #   if (nrow(spell_TSN)>0){
  #     cat("TSNs have some alternatives spellings to review")
  #     #browser()
  #     # for (i in 1:nrow(spell_TSN)){
  #     #   
  #     # }
  #   }
  #   cat(paste0("Alternative spelling review completed in",format(.POSIXct(difftime(Sys.time(), start_spell, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
  # } 
  
  #join results back to original data
  if (doSci) spec_list$SCI_COL_CLN<-NULL
  if (doComm) spec_list$COMM_COL_CLN<-NULL
  
  spec_list_final = spec_list
  spec_list_final$APHIAID_MULTI_FLAG <-FALSE
  spec_list_final$TSN_MULTI_FLAG <-FALSE
  
  if ("APHIAID" %in% codes) {
    
    aphiaids = rbind(definitiveAPHIAID, mysteryAPHIAID)  
    if (nrow(aphiaids[(duplicated(aphiaids$ID, fromLast = FALSE)|duplicated(aphiaids$ID, fromLast = TRUE)),])>0){
      aphiaids_multi = aphiaids[(duplicated(aphiaids$ID, fromLast = FALSE)|duplicated(aphiaids$ID, fromLast = TRUE)),]
    }
    aphiaids = aphiaids[,c("ID","CODE", "CODE_SVC","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")]
    colnames(aphiaids)[colnames(aphiaids) == 'SUGG_SPELLING'] <- 'CODE_SPELLING'
    colnames(aphiaids) <- sub("CODE", "APHIAID", colnames(aphiaids))
    
    if (exists("aphiaids_multi")) aphiaids = aphiaids[!(aphiaids$ID %in% aphiaids_multi$ID),]
    
    spec_list_final = merge(spec_list_final, aphiaids, by="ID", all.x = T)
    if (exists("aphiaids_multi")) {
      spec_list_final[(spec_list_final$ID %in% aphiaids_multi$ID), "APHIAID_MULTI_FLAG"] <-TRUE
    }
  }
  
  if ("TSN" %in% codes) {
    tsns = rbind(definitiveTSN, mysteryTSN)
    if (nrow(tsns[(duplicated(tsns$ID, fromLast = FALSE)|duplicated(tsns$ID, fromLast = TRUE)),])>0){
      tsn_multi = tsns[(duplicated(tsns$ID, fromLast = FALSE)|duplicated(tsns$ID, fromLast = TRUE)),]
    }
    tsns = tsns[,c("ID","CODE", "CODE_SVC","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")]
    colnames(tsns)[colnames(tsns) == 'SUGG_SPELLING'] <- 'CODE_SPELLING'
    colnames(tsns) <- sub("CODE", "TSN", colnames(tsns))
    
    if (exists("tsn_multi")) tsns = tsns[!(tsns$ID %in% tsn_multi$ID),]
    spec_list_final = merge(spec_list_final, tsns, by="ID", all.x = T)
    if (exists("tsn_multi")) {
      spec_list_final[(spec_list_final$ID %in% tsn_multi$ID), "TSN_MULTI_FLAG"] <-TRUE
    }
  }
  
  
  spec_list_final$ID<-NULL
  spec_list_final$ID_SRC <- paste0("Mar.odissupport::getTaxaIDs.R (v",utils::packageDescription('Mar.odissupport')$Version,")") 
  
  if (exists("aphiaids_multi") & exists("tsn_multi")){
    # aphiaids_multiKEEP<<-aphiaids_multi
    # tsn_multiKEEP<<-tsn_multi
    # spec_list_finalKEEP<<-spec_list_final
    multi_final = rbind(aphiaids_multi, tsn_multi)
  } else if (exists("aphiaids_multi")){
    multi_final = aphiaids_multi
  }else if (exists("tsn_multi")){
    multi_final = tsn_multi
  }else{
    multi_final="none"
  }
  if (class(multi_final) == 'data.frame'){
    cat(paste0("Multiple codes were found for some species\n"), file = logName, append = TRUE)
    colnames(multi_final)[colnames(multi_final) == 'CODE'] <- 'CODE_SUGG'
    multi_final = merge(spec_list, multi_final[,-which(colnames(multi_final) %in% c("SCI_COL_CLN","COMM_COL_CLN"))], by="ID", all.y = T)
    multi_final$ID_SRC <- paste0("Mar.odissupport::getTaxaIDs.R (v",utils::packageDescription('Mar.odissupport')$Version,")") 
    
    multi_final$ID<-NULL
  }
  
  res = list(spec_list_final, multi_final)
  cat(paste0("Completed in ",format(.POSIXct(difftime(Sys.time(), start, units="secs"),tz="GMT")),"\n"), file = logName, append = TRUE)
  cat("####################################################\n", file = logName, append = TRUE)
  cat("Done")
  return(res)
}