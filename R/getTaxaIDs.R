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
  
  reviewSpelling <-function(df){
    #this function just drops suggested spelling that match what we already have
    if ("SCI_COL_CLN" %in% names(df) & nrow(df)>0)  df[which(df$SCI_COL_CLN==trimws(gsub(df$SUGG_SPELLING, pattern = "\\(SCIENTIFIC\\)", replacement = ""))),"SUGG_SPELLING"]<-NA
    if ("COMM_COL_CLN" %in% names(df) & nrow(df)>0) df[which(df$COMM_COL_CLN==trimws(gsub(df$SUGG_SPELLING, pattern = "\\(COMMON\\)", replacement = ""))),"SUGG_SPELLING"]<-NA
    return(df)
  }
  
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
  # 
  # if ("TSN" %in% codes) {
  #   #knownAPHIAID<-aphia_uncertain<-data.frame()
  #   knownAPHIAID<-data.frame()
  # }
  # if ("APHIAID" %in% codes) {
  #   mysteryAPHIAID<-definitive
  #   knownTSN<-data.frame()
  #   #knownTSN<-tsn_uncertain<-data.frame()
  # }
    mysteryAPHIAID<-mysteryTSN<-definitive
    knownAPHIAID<-uncertainAPHIAID<-knownTSN<-uncertainTSN<-data.frame()
  rm(definitive)

  if("APHIAID" %in% codes){
    start_AphiaID = Sys.time()
    cat(paste0("AphiaID checks starting\n"), file = logName, append = TRUE)
    if (doSci) {
      #TAXIZE - science
      if (nrow(mysteryAPHIAID)>0){
        tmp_sci_taxize = do_taxize(mysteryAPHIAID,"SCI_COL_CLN",logName,'scientific')
        defCheck = assignDefinitive(df = tmp_sci_taxize, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        if (!exists("definitiveAPHIAID")){
          definitiveAPHIAID<-newdefinitive
        }else{
          definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
      #WORRMS - science
      if (nrow(mysteryAPHIAID)>0){
        tmp_sci_worrms = do_worrms(mysteryAPHIAID,"SCI_COL_CLN",logName,'scientific')
        defCheck = assignDefinitive(df = tmp_sci_worrms, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        if (!exists("definitiveAPHIAID")){
          definitiveAPHIAID<-newdefinitive
        }else{
          definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }

    }
    if (doComm) {
      #TAXIZE - common
      if (nrow(mysteryAPHIAID)>0){
        tmp_comm_taxize = do_taxize(mysteryAPHIAID,"COMM_COL_CLN",logName,'common')
        defCheck = assignDefinitive(df = tmp_comm_taxize, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        if (!exists("definitiveAPHIAID")){
          definitiveAPHIAID<-newdefinitive
        }else{
          definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
      #WORRMS - common

      if (nrow(mysteryAPHIAID)>0){
        tmp_comm_worrms = do_worrms(mysteryAPHIAID,"COMM_COL_CLN",logName,'common')
        defCheck = assignDefinitive(df = tmp_comm_worrms, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        if (!exists("definitiveAPHIAID")){
          definitiveAPHIAID<-newdefinitive
        }else{
          definitiveAPHIAID <- unique(rbind(definitiveAPHIAID[definitiveAPHIAID$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
    }
    cat(paste0("AphiaID checks completed in ",format(.POSIXct(difftime(Sys.time(), start_AphiaID, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
    if(exists("definitiveAPHIAID")) definitiveAPHIAID<-reviewSpelling(definitiveAPHIAID)
    mysteryAPHIAID<-reviewSpelling(mysteryAPHIAID)
    #if a mystery code matches the spelling of the sent value, we'll use  that one and drop the others
    uncertainAPHIAID <- mysteryAPHIAID[is.na(mysteryAPHIAID$SUGG_SPELLING),] 
    mysteryAPHIAID <- mysteryAPHIAID[!(mysteryAPHIAID$ID %in% uncertainAPHIAID$ID),]
    mysteryAPHIAID<-rbind(uncertainAPHIAID,mysteryAPHIAID)
    rm(uncertainAPHIAID)
  }

  if("TSN" %in% codes){
    start_TSN = Sys.time()
    cat(paste0("TSN checks starting\n"), file = logName, append = TRUE)
    if ("APHIAID" %in% codes){
      recs_aphiaid = unique(definitiveAPHIAID[!is.na(definitiveAPHIAID$CODE),])
      aphiaTSNRes =   do_worrmsTSN(recs_aphiaid$CODE,recs_aphiaid, logName=logName)
      #some field renaming below to prevent a warning
      names(aphiaTSNRes)[names(aphiaTSNRes) == 'CODE'] <- 'TSN'
      tsnFields = c(fields, "CODE")
      aphiaTSNRes = merge(recs_aphiaid[,tsnFields],aphiaTSNRes, by.x = "CODE", by.y="APHIAID")
      aphiaTSNRes$CODE<-NULL
      names(aphiaTSNRes)[names(aphiaTSNRes) == 'TSN'] <- 'CODE'
      defCheck = assignDefinitive(df = aphiaTSNRes, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
           newdefinitive= defCheck[[1]]
           mysteryTSNFromAPHIAID= defCheck[[2]]
           if (!exists("definitiveTSN")){
             definitiveTSN<-newdefinitive
           }else{
             definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
           }
           #this ensures that the mystery values include the ones that didn't get sent have an aphiaid
           mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% definitiveTSN$ID),]
           mysteryTSN = mysteryTSN[!(mysteryTSN$ID %in% mysteryTSNFromAPHIAID$ID),]
           mysteryTSN =  rbind(mysteryTSN, mysteryTSNFromAPHIAID)

    }

    if (doSci) {
      # ritis - science
      if (nrow(mysteryTSN)>0){
        tmp_sci_ritis = do_ritis(mysteryTSN,"SCI_COL_CLN",logName,'scientific')
        defCheck = assignDefinitive(df = tmp_sci_ritis, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryTSN= defCheck[[2]]
        if (!exists("definitiveTSN")){
          definitiveTSN<-newdefinitive
        }else{
          definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
    }
     if (doComm) {
      #ritis - science
      if (nrow(mysteryTSN)>0){
        tmp_comm_ritis = do_ritis(mysteryTSN,"COMM_COL_CLN",logName,'common')
        defCheck = assignDefinitive(df = tmp_comm_ritis, masterList = spec_list[,names(spec_list) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")])
        newdefinitive= defCheck[[1]]
        mysteryTSN= defCheck[[2]]
        if (!exists("definitiveTSN")){
          definitiveTSN<-newdefinitive
        }else{
          definitiveTSN <- unique(rbind(definitiveTSN[definitiveTSN$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
    }
    cat(paste0("TSN checks completed in ",format(.POSIXct(difftime(Sys.time(), start_TSN, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
    if(exists("definitiveTSN")) definitiveTSN<-reviewSpelling(definitiveTSN)
    mysteryTSN<-reviewSpelling(mysteryTSN)
    #if a mystery code matches the spelling of the sent value, we'll use  that one and drop the others
    uncertainTSN <- mysteryTSN[is.na(mysteryTSN$SUGG_SPELLING),] 
    mysteryTSN <- mysteryTSN[!(mysteryTSN$ID %in% uncertainTSN$ID),]
    mysteryTSN<-rbind(uncertainTSN,mysteryTSN)
    rm(uncertainTSN)
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