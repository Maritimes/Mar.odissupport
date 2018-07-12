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
#' @family speciesCodes
#' @importFrom stats setNames
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
getTaxaIDs <- function(spec_list = NULL,
                       sci_col = NULL,
                       comm_col = NULL,
                       sci_Filts = NULL,
                       comm_Filts = NULL,
                       codes = c("APHIAID", "TSN")) {
  log_nm = "getTaxaIDs.log"
  log_con <- file(log_nm)
  cat(paste0("Initiated at ",format(Sys.time(), "%Y-%m-%d %H:%M"),"\n"), file = log_con) 
  
  cat(paste0("A log file has been generated at ",file.path(getwd(),log_nm),".\nIn case of failure, that file will reveal the last successful operation"))
  #filters are the same, but might be handy to be able to remove things from
  #common names but not scientific, and vice versa
  allFilts <- c("BAIT", "DIGESTED","UNIDENTIFIED PER","UNIDENTIFIED SPECIES",
                "UNID (FISH|REMAINS)+", "REMAINS","SURVEY","FSRS -","RESERVED",
                "PURSE","^FISH( AND| \\,|$)","\\,?\\s?EGG(S?)-?","\\s?LARVAE",
                "INVERTEBRATE","WATER","FLUID","^SAND$","EGGS",
                "INORGANIC DEBRIS","MIXED","MUCUS","OPERCULUM","^SHARK$")
  
  commFilts <- c(comm_Filts,"([^']\\b[SP]{1,3}\\.?$)",
                 "([^']\\b[a-zA-Z]{1,2}\\.?$)","SHARK ","^FISH$","/")
  
  sciFilts <- c(sci_Filts, "WHALE","CETACEAN","/","CRAB", "LOBSTER","SHRIMP",
                "IRISH MOSS","SHARK","COD WORM","SEA CORALS","SKATE","OBSOLETE",
                "FINFISHES","GROUNDFISH","PELAGIC FISH","\\bAND\\b","SAND TUBE",
                "UNIDENTIFIED",",")
  spec_list$ID <- seq.int(nrow(spec_list))
  
  spec_list$SCI_COL_CLN = NA
  doSci = F
  spec_list$COMM_COL_CLN = NA
  doComm = F

  if (!is.null(sci_col)) {
    #remove whitespace
    spec_list$SCI_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list[, sci_col]))
    #remove recs matching allFilts and sci_filts
    spec_list[grepl(x=spec_list[, sci_col],ignore.case = TRUE, pattern = paste(c(allFilts, sciFilts), collapse = "|")),"SCI_COL_CLN"]<-NA
    
    #remove bad bits, but retain the rest of the string)

    #SCI names never have brackets - get rid of them, and everything they contain
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\(.*?\\)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\(.*?\\)"),"SCI_COL_CLN"],
           pattern = "\\(.*?\\)",replacement = "") 
    
    #(NS) or NS
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"SCI_COL_CLN"],
           pattern = "(\\(NS\\)|\\bNS)",replacement = "") 
    #SP, SP., SPP and SPP.
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "([^']\\b[SP]{1,3}\\.?$)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "([^']\\b[SP]{1,3}\\.?$)"),"SCI_COL_CLN"],
           pattern = "([^']\\b[SP]{1,3}\\.?$)",replacement = "") 
    
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\b[a-zA-Z]{1,2}\\."),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\b[a-zA-Z]{1,2}\\."),"SCI_COL_CLN"],
           pattern = "\\b[a-zA-Z]{1,2}\\.",replacement = "") 
    
    #final remove whitespace and drop really short records
    spec_list[,"SCI_COL_CLN"]<-  gsub("(^\\s+)|(\\s+$)", "", spec_list$SCI_COL_CLN)
    spec_list[which(nchar(spec_list$SCI_COL_CLN)<4),"SCI_COL_CLN"]<-NA
    doSci = T
  } 
  if (!is.null(comm_col)) {
    #remove whitespace
    spec_list$COMM_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list[, comm_col]))
    #remove recs matching allFilts and sci_filts
    spec_list[grepl(x=spec_list[, comm_col],ignore.case = TRUE, pattern = paste(c(allFilts, commFilts), collapse = "|")),"COMM_COL_CLN"]<-NA
    
    #remove bad bits, but retain the rest of the string)
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(-|,|\\s)?UNIDENTIFIED.*"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(-|,|\\s)?UNIDENTIFIED.*"),"COMM_COL_CLN"],
           pattern = "(-|,|\\s)?UNIDENTIFIED.*",replacement = "") 
     spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*"),"COMM_COL_CLN"]<-
       gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*"),"COMM_COL_CLN"],
            pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*",replacement = "")
    #(NS) or NS
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"COMM_COL_CLN"],
           pattern = "(\\(NS\\)|\\bNS)",replacement = "") 
    
    #stupid "FISH" record
    spec_list[grepl(x = spec_list[,comm_col],ignore.case = T,pattern = "UNID\\. FISH"),"COMM_COL_CLN"]<-NA
    #final remove whitespace and drop really short records
    spec_list$COMM_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list$COMM_COL_CLN))
    spec_list[which(nchar(spec_list$COMM_COL_CLN)<4),"COMM_COL_CLN"]<-NA
    doComm = T
  } 
  
  cat(paste0("Input data successfully filtered\n"), file = "getTaxaIDs.log", append = TRUE)
  definitive = spec_list[, c("ID", "SCI_COL_CLN", "COMM_COL_CLN")]
  cols = c("CODE","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","CODE_SRC","SUGG_SPELLING")
  definitive = cbind(definitive, setNames( lapply(cols, function(x) x=NA), cols) )
  
  definitiveTSN = definitive
  mysteryTSN = definitive
  
  definitiveAPHIAID = definitive
  mysteryAPHIAID = definitive
  rm(definitive)
  
  aphia_known = data.frame()
  aphia_uncertain = data.frame()
  aphia_mystery = data.frame()
  tsn_known = data.frame()
  tsn_uncertain = data.frame()
  tsn_mystery = data.frame()
  
  if("APHIAID" %in% codes){
    
    cat(paste0("AphiaID checks starting\n"), file = "getTaxaIDs.log", append = TRUE)
    theAphiaIDs <- getAphiaIDs(mysteryAPHIAID = mysteryAPHIAID, 
                               masterList = spec_list[,c("ID", "SCI_COL_CLN", "COMM_COL_CLN"),], 
                               doSci=doSci, doComm=doComm)
    aphia_known = theAphiaIDs[[1]]
    
    aphia_uncertain = theAphiaIDs[[2]][!is.na(theAphiaIDs[[2]]$CODE),]
    #remove unnecessary spelling tips
    #sci
    if (nrow(aphia_known[which(toupper(trimws(aphia_known$SCI_COL_CLN))==toupper(trimws(aphia_known$SUGG_SPELLING))),])>0)
      aphia_known[which(toupper(trimws(aphia_known$SCI_COL_CLN))==toupper(trimws(aphia_known$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    if (nrow(aphia_known[which(toupper(trimws(aphia_known$COMM_COL_CLN))==toupper(trimws(aphia_known$SUGG_SPELLING))),])>0)
      aphia_known[which(toupper(trimws(aphia_known$COMM_COL_CLN))==toupper(trimws(aphia_known$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    #common
    if (nrow(aphia_uncertain[which(toupper(trimws(aphia_uncertain$SCI_COL_CLN))==toupper(trimws(aphia_uncertain$SUGG_SPELLING))),])>0)
      aphia_uncertain[which(toupper(trimws(aphia_uncertain$SCI_COL_CLN))==toupper(trimws(aphia_uncertain$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    if (nrow(aphia_uncertain[which(toupper(trimws(aphia_uncertain$COMM_COL_CLN))==toupper(trimws(aphia_uncertain$SUGG_SPELLING))),])>0)
      aphia_uncertain[which(toupper(trimws(aphia_uncertain$COMM_COL_CLN))==toupper(trimws(aphia_uncertain$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    
    aphia_mystery = theAphiaIDs[[2]][is.na(theAphiaIDs[[2]]$CODE),]
    if (nrow(aphia_mystery)>0) 
      aphia_mystery[is.na(aphia_mystery$CODE),c("CODE_SVC","CODE_TYPE","CODE_SRC","SUGG_SPELLING")]<-NA
  }
  
  if("TSN" %in% codes){
    
    cat(paste0("TSN checks starting\n"), file = "getTaxaIDs.log", append = TRUE)
    if ("APHIAID" %in% codes) {
      knownAphias <- theAphiaIDs[[1]][!is.na(theAphiaIDs[[1]]$CODE),]
      if (nrow(knownAphias)<1)knownAphias <-NULL
    }else{
      knownAphias <-NULL
    }
    theTSNs <- getTSNs(mysteryTSN = mysteryTSN, 
                       masterList = spec_list[,c("ID", "SCI_COL_CLN", "COMM_COL_CLN")], 
                       knownAphias = knownAphias, 
                       doSci=doSci, doComm=doComm)
    tsn_known = theTSNs[[1]]
    tsn_uncertain = theTSNs[[2]][!is.na(theTSNs[[2]]$CODE),]
    #remove unnecessary spelling tips
    if (nrow(tsn_known[which(toupper(trimws(tsn_known$SCI_COL_CLN))==toupper(trimws(tsn_known$SUGG_SPELLING))),])>0)
      tsn_known[which(toupper(trimws(tsn_known$SCI_COL_CLN))==toupper(trimws(tsn_known$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    if (nrow(tsn_known[which(toupper(trimws(tsn_known$COMM_COL_CLN))==toupper(trimws(tsn_known$SUGG_SPELLING))),])>0)
      tsn_known[which(toupper(trimws(tsn_known$COMM_COL_CLN))==toupper(trimws(tsn_known$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    if (nrow(tsn_uncertain[which(toupper(trimws(tsn_uncertain$SCI_COL_CLN))==toupper(trimws(tsn_uncertain$SUGG_SPELLING))),])>0)
      tsn_uncertain[which(toupper(trimws(tsn_uncertain$SCI_COL_CLN))==toupper(trimws(tsn_uncertain$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    if (nrow(tsn_uncertain[which(toupper(trimws(tsn_uncertain$COMM_COL_CLN))==toupper(trimws(tsn_uncertain$SUGG_SPELLING))),])>0)
      tsn_uncertain[which(toupper(trimws(tsn_uncertain$COMM_COL_CLN))==toupper(trimws(tsn_uncertain$SUGG_SPELLING))),"SUGG_SPELLING"]<-NA
    
    tsn_mystery = theTSNs[[2]][is.na(theTSNs[[2]]$CODE),]
    if (nrow(tsn_mystery)>0) 
      tsn_mystery[,c("CODE_SVC","CODE_TYPE","CODE_SRC","SUGG_SPELLING")]<-NA
  }
  
  if ("TSN" %in% codes & "APHIAID" %in% codes ){
    #this component tries to find aphiaids using any suggested spellings that 
    #came from finding the TSN
    dblCheck = rbind(aphia_uncertain,aphia_mystery)
    if(nrow(dblCheck)>0){
      dblCheck = dblCheck[dblCheck$ID %in% tsn_known[!is.na(tsn_known$SUGG_SPELLING) ,"ID"],]
      if(nrow(dblCheck)>0){
        correctors = tsn_known[tsn_known$ID %in% dblCheck$ID,c("ID","SUGG_SPELLING","CODE_SVC")]
        if (nrow(correctors)>0){
          cat(paste0("AphiaID re-check initiated (using suggested spellings)\n"), file = "getTaxaIDs.log", append = TRUE)
          correctors$SCI_COL_CLN<-correctors$COMM_COL_CLN<-correctors$SUGG_SPELLING
          correctors$SUGG_SPELLING<-NULL
          dblCheckTAXSCI =   chk_taxize(correctors, "SCI_COL_CLN",searchtype = 'scientific')
          dblCheckWORRMSSCI =   chk_worrms(correctors, "SCI_COL_CLN",searchtype = 'scientific')
          newres = rbind(dblCheckTAXSCI,dblCheckWORRMSSCI) 
          newres=newres[newres$CODE_DEFINITIVE %in% TRUE,]
          newres = merge(newres, correctors[,c("ID","CODE_SVC")], by = "ID")
          newres[,"CODE_SRC"]<-paste0("scientific via (",newres$CODE_SVC.y,")")
          newres$CODE_SVC.y<-NULL
          colnames(newres)[colnames(newres) == 'CODE_SVC.x'] <- 'CODE_SVC'
          defCheck = assignDefinitive(df = newres, masterList = spec_list[,c("ID", "SCI_COL_CLN", "COMM_COL_CLN"),])
          newdefinitive= defCheck[[1]]
          newAphia_mystery= defCheck[[2]]
          if (!exists("aphia_known")){
            aphia_known<-newdefinitive
          }else{
            aphia_known <- unique(rbind(aphia_known,newdefinitive))
          }
          aphia_mystery = rbind(aphia_mystery,newAphia_mystery)
          aphia_mystery = aphia_mystery[!aphia_mystery$ID %in% aphia_known$ID,]
          rm(defCheck)
          rm(newdefinitive)
          rm(newAphia_mystery)
        }
      }
    }
    
  }
  
  #join results back to original data
  spec_list$SCI_COL_CLN<-NULL
  #spec_list$COMM_COL_CLN<-NULL
  
  spec_list_final = spec_list
  spec_list_final$APHIAID_MULTI_FLAG <-FALSE
  spec_list_final$TSN_MULTI_FLAG <-FALSE
  # aphiaids_multi <- "none"
  # tsns_multi <- "none"
  
  if ("APHIAID" %in% codes) {
    aphiaids = rbind(aphia_known, aphia_uncertain, aphia_mystery)  
    if (nrow(aphiaids[(duplicated(aphiaids$ID, fromLast = FALSE)|duplicated(aphiaids$ID, fromLast = TRUE)),])>0){
      aphiaids_multi = aphiaids[(duplicated(aphiaids$ID, fromLast = FALSE)|duplicated(aphiaids$ID, fromLast = TRUE)),]
      aphiaids_multi$CODE_ID<-"APHIAID"
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
    tsns = rbind(tsn_known, tsn_uncertain, tsn_mystery)
    if (nrow(tsns[(duplicated(tsns$ID, fromLast = FALSE)|duplicated(tsns$ID, fromLast = TRUE)),])>0){
      tsn_multi = tsns[(duplicated(tsns$ID, fromLast = FALSE)|duplicated(tsns$ID, fromLast = TRUE)),]
      tsn_multi$CODE_ID<-"TSN"
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
    multi_final = cbind(aphiaids_multi, tsn_multi)
  } else if (exists("aphiaids_multi")){
    multi_final = aphiaids_multi
  }else if (exists("tsn_multi")){
    multi_final = tsn_multi
  }else{
    multi_final="none"
  }
  if (class(multi_final) == 'data.frame'){
    cat(paste0("Multiple codes were found for some species\n"), file = "getTaxaIDs.log", append = TRUE)
    colnames(multi_final)[colnames(multi_final) == 'CODE'] <- 'CODE_SUGG'
    multi_final = merge(spec_list, multi_final, by="ID", all.y = T)
    multi_final$ID<-NULL
    multi_final$SCI_COL_CLN<-NULL
    #multi_final$COMM_COL_CLN<-NULL
  }
  
  res = list(spec_list_final, multi_final)
  
  cat("####################################################\n", file = "getTaxaIDs.log", append = TRUE)
  return(res)
}