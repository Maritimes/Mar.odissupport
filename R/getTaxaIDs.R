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
#'   \item Aphiaid
#'   \item Aphiaid_src - which service(s) provided the Aphiaid value
#'   \item Aphiaid_data - which service(s) provided the Aphiaid value
#'   \item Aphiaid_definitive* - indicates whether or not we are confident about 
#'   AphiaID
#'   \item TSN 
#'   \item TSN_src - which service(s) provided the TSN value
#'   \item TSN_data - which service(s) provided the TSN value
#'   \item TSN_definitive* - indicates whether or not we are confident about 
#'   TSN
#'   }
#' Most webservices can return multiple IDs for a given search term, and the 
#' results are often identified as "accepted" or not. For this function to 
#' identify an ID as being definitive=TRUE, one of the following 3 situations 
#' occurred: 
#' \itemize{
#' \item Only a single result was returned from a service
#' \item Only a single item within a resultset was marked as "Accepted"
#' \item Two independent services both returned the same ID (These may be 
#' considered less definitive, and can be recognized by their <>_src fields 
#' indication multiple sources - e.g. "Both worrms & ritis")
#' }

#' @param spec_list the dataframe containing information to be decoded to TSN and
#' aphiaIDs
#' @param sci_col the name of the column of the dataframe containing the
#' scientific names
#' @param comm_col the name of the column of the dataframe containing the
#' common names
#' @param filterStrings a vector of regex values that you might want to clean out
#' of your scientific and common names prior to sending them to a web service.  
#' For example, some Maritimes names inlude "(NS)", which will prevent services
#' from finding matches, By adding "\\(NS\\)" (escaping the brackets and periods 
#' with slashes), we ensure that the results will be as clean as possible prior 
#' to searching. By default, the following values will be filtered:
#' \itemize{
#' \item \\(NS\\)
#' \item SP\\.
#' \item S\\.O\\.
#' \item F\\.
#' \item UNIDENTIFIED
#' \item EGGS
#' \item UNID\\.
#' \item \\(ORDER\\)
#' \item \\(SUBORDER\\)
#' }
#' @importFrom utils head
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom ritis search_scientific
#' @importFrom ritis search_common
#' @importFrom taxize itis_acceptname
#' @importFrom taxize get_wormsid
#' @importFrom jsonlite fromJSON
#' @family qc
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
#' @note
getTaxaIDs <- function(spec_list=NULL, sci_col=NULL, comm_col=NULL, filterStrings=NULL){
  filterStrings <- c(filterStrings, "\\(NS\\)", " SP\\.", "S\\.O\\.", "F\\.","UNIDENTIFIED", "EGGS","UNID\\.","\\(ORDER\\)","\\(SUBORDER\\)")
  
  # Initial data preparation ------------------------------------------------
  spec_list_orig = spec_list
  if (is.null(comm_col)){
    spec_list<-data.frame("tmpzzz" = spec_list[,sci_col])
    names(spec_list)[names(spec_list)=="tmpzzz"]<-sci_col
  } else {
    spec_list<-data.frame(spec_list[,c(sci_col, comm_col)])
  }
  
  #' remove potential confusion from CASE sensitivity
  if (nrow(spec_list)==1){
    spec_list = as.data.frame(t(as.data.frame(sapply(spec_list, toupper))))
  }else{
    spec_list = as.data.frame(sapply(spec_list, toupper))
  }
  resultFormat = data.frame(internal_name = character(), suggested_name = character(), search_type = character(),src= character(), ID= integer())
  
  #'_src is where we got the ID (e.g. worrms)
  #'_data is what we used to get the id (e.g. scientific/common name, or another ID)
  #multi_Codes will hold alternative id values for species where we found multiple matches
  multi_Codes <<- data.frame(suggested_name = character(), TSN= integer(), TSN_internal_name = character(), TSN_src= character(), TSN_data = character(),
                             AphiaID_internal_name = character(), AphiaID_src= character(), AphiaID_data = character(), AphiaID= integer() )
  mysterySpec <-spec_list
  #vector of string we can try dropping
  
  definitive <- spec_list
  names(definitive)[names(definitive)==sci_col]<-"SPEC"
  names(definitive)[names(definitive)==comm_col]<-"COMM"
  definitive$AphiaID_Definitive <- NA
  definitive$AphiaID <- NA
  definitive$AphiaID_src <- NA
  definitive$AphiaID_data <- NA
  definitive$TSN_Definitive <- NA
  definitive$TSN <- NA
  definitive$TSN_src <- NA
  definitive$TSN_data <- NA
  df_res<-spec_list
  # Basic helper functions --------------------------------------------------
  
  sapply_pb <- function(X, FUN, ...)
  {
    #'this is a progress bar function I can wrap around stuff
    #'think I stole it from https://ryouready.wordpress.com/2010/01/11/
    #'progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
    }
    res <- sapply(X, wrapper, ...)
    close(pb)
    res
  }
  
  multiCheck<-function(multi_Codes){
    #update multi_codes to hold only those species for which we don't have a definitve value 
    theseMulti = multi_Codes[(multi_Codes$SPEC %in% definitive[definitive$Is_Definitive==FALSE,"SPEC"]),]
    assign("multi_Codes", theseMulti, envir = .GlobalEnv)
    return(theseMulti)
  }
 
  # Functions for checking web services -------------------------------------
  chk_taxize <- function(filterStrings, searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use taxize functions to check stuff
    if (searchtype == 'scientific'){
      colname = "SPEC"
    }else{
      colname = "COMM"
    }
    #vector of string we can try dropping
    if (grepl(paste(filterStrings,collapse="|"),searchterm)==TRUE){
      tmpTerm = trimws(gsub(paste(filterStrings,collapse="|"),  "", searchterm))
    }else{
      tmpTerm = searchterm
    }
    this <- tryCatch(
      {
        taxize::get_wormsid(query=tmpTerm, searchtype=searchtype, ask=ask, verbose=verbose, accepted = TRUE,rows = 1, messages = F)
      },
      error=function(cond){
      }
    )
    if (is.null(this)){
      cmt = "not found"
      ID = NA
    }else {
      ID =as.vector(this)
      if (attr(this,"match")=="found" & attr(this,"multiple_matches") == TRUE) {
        cmt = "multi match - kept 1st accepted"
        mult = spec_list[spec_list[,colname] == searchterm,]
        mult$AphiaID = ID
        mult$AphiaID_src = "Taxize"
        mult$AphiaID_data = searchtype
        mult$AphiaID_internal_name = searchterm
        mult$suggested_name = NA #can't get names back
        mult$TSN_internal_name = NA
        mult$TSN_src = NA
        mult$TSN_data = NA
        mult$TSN = NA
        multi_Codes <- rbind(.GlobalEnv$multi_Codes, mult)
        assign("multi_Codes", multi_Codes, envir = .GlobalEnv)
        
        #can't capture all other matches - only can get first
      } else if (attr(this,"match")=="found" & attr(this,"multiple_matches") == FALSE) {
        cmt = "multi match - kept only accepted"
      }else{
        cmt = "multi match - no ids returned"
      }
    }
    this = as.data.frame(cbind("SEARCHTERM" = searchterm, "ID_taxize"=ID,"MATCH_taxize"= cmt))
    if (NROW(this[grep("NA due to ask=FALSE", this$MATCH_WORMS),])>0)  this[grep("NA due to ask=FALSE", this$MATCH_WORMS),]$MATCH_WORMS<-"multi match"
    return(this)
  }
  
  chk_worrms <- function(filterStrings, searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    if (searchtype == 'scientific'){
      colname = "SPEC"
    }else{
      colname = "COMM"
    }
    if (grepl(paste(filterStrings,collapse="|"),searchterm)==TRUE){
      tmpTerm = trimws(gsub(paste(filterStrings,collapse="|"),  "", searchterm))
    }else{
      tmpTerm = searchterm
    }
    resName = "scientificname"
    if (searchtype == 'scientific'){
      thisi <- tryCatch(
        {
          worrms::wm_records_name(tmpTerm, marine_only = T)
        },
        error=function(cond){
        }
      )
    }else{
      thisi <- tryCatch(
        {
          worrms::wm_records_common(tmpTerm, marine_only = T)
        },
        error=function(cond){
        }
      )
    }
    crapRec = data.frame("SEARCHTERM" = searchterm, "ID_worrms" = NA, "MATCH_worrms" = "not found")
    if (is.null(thisi)) {
      thisi <- crapRec
    } else if (nrow(thisi)==0) {
      thisi <- crapRec
    }else if (nrow(thisi)==1) {
      thisi = cbind("SEARCHTERM"= searchterm,"ID_worrms" = thisi$AphiaID, "MATCH_worrms" = "single match")
    }else{
      if (nrow(thisi[thisi$status=="accepted",])==1){
        #maybe need merge thisi with speclist?
        thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi[thisi$status=="accepted","AphiaID"], "MATCH_worrms" = "multi match - kept only accepted")
        colnames(thisi)[colnames(thisi)=="AphiaID"] <- "ID_worrms"
      } else {
        mult = thisi[,c(resName,"AphiaID")]
        colnames(mult)<-c("suggested_name","AphiaID")
        mult$AphiaID_src = "worrms"
        mult$joiner = searchterm
        mult$AphiaID_internal_name = searchterm
        mult$AphiaID_data = searchtype
        mult$TSN_internal_name = NA
        mult$TSN_src = NA
        mult$TSN_data = NA
        mult$TSN = NA
        spec_list$joinerx<-spec_list[,colname]
        mult = merge(spec_list, mult, by.x="joinerx", by.y="joiner")
        mult$joinerx<-NULL
        mult$joiner = NULL
        multi_Codes <- rbind(.GlobalEnv$multi_Codes, mult)
        assign("multi_Codes", multi_Codes, envir = .GlobalEnv)
        if (nrow(thisi[thisi$status=="accepted",])>1) {
          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - kept 1st accepted")
        }else {
          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - none accepted; retained first")
        }
      }
    }
    #res=list(thisi,multi_Codes)
    return(thisi)
  }
  
  chk_ritis <- function(filterStrings, searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use ritis functions to check stuff
    this = resultFormat
    if (grepl(paste(filterStrings,collapse="|"),searchterm)==TRUE){
      tmpTerm = trimws(gsub(paste(filterStrings,collapse="|"),  "", searchterm))
    }else{
      tmpTerm = searchterm
    }
    if (searchtype == 'scientific'){
      resName = "combinedName"
      colname = "SPEC"
      thisi =  as.data.frame(ritis::search_scientific(tmpTerm, wt = "json", raw = FALSE))
    }else{
      resName = "commonName"
      colname = "COMM"
      thisi = as.data.frame(ritis::search_common(tmpTerm, wt = "json", raw = FALSE))
    }
    if (nrow(thisi)==0) {
      thisi=data.frame("SEARCHTERM" = searchterm, "ID_ritis" = NA, "MATCH_ritis" = "not found")
    }else if (nrow(thisi)==1) {
      thisi = cbind("SEARCHTERM"= searchterm,"ID_ritis" = thisi$tsn, "MATCH_ritis" = "single match")
    }else{
      thisi[,resName]<-toupper(thisi[,resName])
      mult = thisi[,c(resName,"tsn")]
      colnames(mult)<-c("suggested_name","TSN")
      mult$TSN_internal_name = searchterm
      mult$joiner = searchterm
      mult$TSN_src = "ritis"
      mult$TSN_data = searchtype
      mult$AphiaID_internal_name = NA
      mult$AphiaID_src = NA
      mult$AphiaID_data = NA
      mult$AphiaID = NA
      mult = merge(mult,spec_list, by.x = "joiner", by.y = colname) 
      mult$TSN<-mult$TSN.x
      mult$TSN.x<-NULL
      mult$TSN.y<-NULL
      mult$joiner = NULL
      #multi_Codes <<- rbind(multi_Codes, mult)
      if (length(thisi[thisi[,resName] == searchterm,"tsn"])==1){
        thisi = cbind("SEARCHTERM"=searchterm, "ID_ritis" = thisi[thisi[,resName] == searchterm,"tsn"], "MATCH_ritis" ="multi match - retained identical spelling")
      }else{
        thisi = cbind("SEARCHTERM"=searchterm, "ID_ritis" =thisi$tsn[1], "MATCH_ritis" ="multi match - no identical spelling; retained first")
      }
    }
    return(thisi)
  }
  
  chk_worrmsTSN<-function(AphiaIDs){
    thisi <- tryCatch(
      {
        worrms::wm_external_(id = as.integer(AphiaIDs))
      },
      error=function(cond){
      }
    )
    return(thisi)
  }
  
  # Main Functions ----------------------------------------------------------
  assignDefinitive<-function(curDefinitive, results = NULL, searchtype = NULL){
    if (searchtype =="scientific"){
      searchtypeVal = "scientific name"
      colname = "SPEC"
    }else{
      searchtypeVal = "common name"
      colname = "COMM"
    }
    if(any(grep("taxize",colnames(results))) |(any(grep("worrms",colnames(results))) )){
      theID = "AphiaID"
      theIDComment = "AphiaID_src"
      theIDUsed = "AphiaID_data"
      #defTableName = "definitive_Aphia"
      isDefFlag = "AphiaID_Definitive"
      
      otherTheID = "TSN"
      otherTheIDComment = "TSN_src"
      otherTheIDUsed = "TSN_data"
      #otherDefTableName = "definitive_TSN"
      otherIsDefFlag = "TSN_Definitive"
      if (any(grep("taxize",colnames(results)))){
        results_ID = "ID_taxize"
        results_MATCH = "MATCH_taxize"
        results_type = "Taxize"
      }else if (any(grep("worrms",colnames(results)))){
        results_ID = "ID_worrms"
        results_MATCH = "MATCH_worrms"
        results_type = "Worrms"
      }
    }else if (any(grep("ritis",colnames(results)))){
      theID = "TSN"
      theIDComment = "TSN_src"
      theIDUsed = "TSN_data"
      #defTableName = "definitive_TSN"
      isDefFlag = "TSN_Definitive"
      
      otherTheID = "AphiaID"
      otherTheIDComment = "AphiaID_src"
      otherTheIDUsed = "AphiaID_data"
      #otherDefTableName = "definitive_Aphia"
      otherIsDefFlag = "AphiaID_Definitive"
      
      results_ID = "ID_ritis"
      results_MATCH = "MATCH_ritis"
      results_type = "ritis"
    }
    excludes <- c(theID,theIDComment,theIDUsed,isDefFlag)
    keeperStrings = c("single match","only accepted")
    results[theID] = NA
    results[theIDUsed] = NA
    results[theIDComment] = NA
    results[isDefFlag]=NA
    if (nrow(results[!is.na(results[,results_ID]),])>0){
      #for those which we have results - carry over the ID, the src and how it was found
      results[!is.na(results[,results_ID]),theIDUsed]<-searchtypeVal
      results[!is.na(results[,results_ID]),theIDComment]<-results_type
      results[!is.na(results[,results_ID]) & grepl(paste(keeperStrings,collapse="|"),results[,results_MATCH]),isDefFlag]<-TRUE
      results[!is.na(results[,results_ID]) & !(results[,isDefFlag] %in% TRUE),isDefFlag]   <-FALSE
      results[!is.na(results[,results_ID]),][theID] <- results[!is.na(results[,results_ID]),][results_ID]
    }
    
    results <-results[,c("SEARCHTERM",theID,theIDComment,theIDUsed, isDefFlag)]
    names(results)[names(results) == 'SEARCHTERM'] <- colname
    
    multimatches = 
    
    replacers = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(FALSE)),][c(colname,theID)], 
                      results[results[,isDefFlag] %in% TRUE,], by=colname)
    replacers[,'thisID'] <- replacers[paste0(theID,".y")]
    doubleblind = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(FALSE)),][c(colname,theID)], 
                        results[results[,isDefFlag] %in% FALSE,], by=colname)
    doubleblind[,'thisID'] <- doubleblind[paste0(theID,".y")]
    newFind = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(NA)),][c(colname,theID)], 
                    results[results[,isDefFlag] %in% c(TRUE,FALSE),], by=colname)
    newFind[,'thisID'] <- newFind[paste0(theID,".y")]
    
      #if the new codes are definitive, we will replace the existing (not-definitive) ones entirely 
      if (nrow(replacers)>0){
        curDefinitive[curDefinitive[,colname] %in% replacers[,colname],c(theIDComment,theIDUsed, isDefFlag)]<-
          replacers[c(theIDComment,theIDUsed,isDefFlag)]
      }
      #if we have already have a tentative code, and our new result set gives the same code,
      #we will call it definitive        
       if (nrow(doubleblind)>0){
         for (i in 1:nrow(doubleblind)){
                           curDefinitive[curDefinitive[,isDefFlag] %in% FALSE & curDefinitive[,colname] %in% doubleblind[,colname],][i,theIDComment]<-
              paste("Both",curDefinitive[curDefinitive[,isDefFlag] %in% FALSE & curDefinitive[,colname] %in% doubleblind[,colname],theIDComment][i], "&",doubleblind[i,theIDComment])  #c(theIDComment,theIDUsed)]<-NA
         }
         curDefinitive[curDefinitive[,isDefFlag] %in% FALSE & curDefinitive[,colname] %in% doubleblind[,colname],isDefFlag]<-TRUE
       } 
      #just adding new stuff - we don't have anything yet
      if (nrow(newFind)>0){
        curDefinitive[curDefinitive[,isDefFlag] %in% NA & curDefinitive[,colname] %in% newFind[,colname],][,c(theID,theIDComment,theIDUsed, isDefFlag)]<-
          newFind[,c("thisID",theIDComment,theIDUsed, isDefFlag)]
      }
    return(curDefinitive)
  }   
  checkAndCompareTSN<-function(filterStrings=filterStrings, definitive=definitive, df=NULL, searchtype = NULL){
    if (searchtype =="scientific"){
      searchtypeVal = "scientific name"
      colname = "SPEC"
    }else{
      searchtypeVal = "common name"
      colname = "COMM"
    }
    print(paste("...Checking ritis using ",searchtype," names"))
    no_TSN = definitive[definitive$TSN_Definitive %in% c(NA, FALSE),]
    if (nrow(no_TSN)>0) {
      #print(colnames(no_TSN))
      no_TSN = as.data.frame(t(sapply_pb(no_TSN[,colname],chk_ritis, filterStrings = filterStrings, searchtype=searchtype)))
      names(no_TSN) = c("SEARCHTERM","ID_ritis","MATCH_ritis")
      definitive = assignDefinitive(definitive,no_TSN, searchtype)
    }
    return(definitive)
  }
  
  checkAndCompareAphiaID<-function(definitive, df= NULL, searchtype = NULL){
    if (searchtype =="scientific"){
      colname = "SPEC"
    }else{
      colname = "COMM"
    }
    df=merge(df,spec_list)

    no_aphiaID = df[!(df$AphiaID_Definitive %in% TRUE),]
    
    if (nrow(no_aphiaID)>0) {
      print(paste0("...against taxize using ",searchtype," names"))
      no_aphiaID = as.data.frame(t(sapply_pb(no_aphiaID[,colname],chk_taxize, filterStrings = filterStrings, searchtype=searchtype)))
    }
    if (nrow(no_aphiaID)>0) {
      no_aphiaID = data.frame("SEARCHTERM" = toupper(unlist(no_aphiaID$SEARCHTERM)),
                              "ID_taxize" = unlist(no_aphiaID$ID_taxize),
                              "MATCH_taxize" = unlist(no_aphiaID$MATCH_taxize))
      definitive = assignDefinitive(definitive,no_aphiaID, searchtype)
      not_definitive = definitive[!(definitive$AphiaID_Definitive %in% TRUE),]
      if (nrow(not_definitive)>0){
        print(paste0("...against worrms using ",searchtype," names"))
        not_definitive = as.data.frame(t(sapply_pb(not_definitive[,colname],chk_worrms, filterStrings = filterStrings, searchtype=searchtype)))
        names(not_definitive)<-c("SEARCHTERM", "ID_worrms", "MATCH_worrms")
        not_definitive = data.frame("SEARCHTERM" = toupper(unlist(not_definitive$SEARCHTERM)),
                                    "ID_worrms" = unlist(not_definitive$ID_worrms),
                                    "MATCH_worrms" = unlist(not_definitive$MATCH_worrms))
        definitive = assignDefinitive(definitive,not_definitive, searchtype)
      }
    }
    return(definitive)
  }  
  
  
  # Start of the actual processing ------------------------------------------
  print("Looking for AphiaIDs")
  definitive = checkAndCompareAphiaID(definitive, df = definitive,searchtype="scientific")
  missing = definitive[!(definitive$AphiaID_Definitive %in% TRUE),]
  multi_Codes = multiCheck(multi_Codes)
  if (nrow(missing)>0){
    definitive = checkAndCompareAphiaID(definitive, df = missing,searchtype="common")
    multi_Codes = multiCheck(multi_Codes)
  }
  # #have lots of aphias - try to use them to get TSNs
  print("Looking for TSNs")
  print("...Checking worrms using AphiaIDs")
  if (!is.null(definitive)){
    TSNsFromAphiasDef = sapply_pb(definitive$AphiaID,chk_worrmsTSN)
    TSNsFromAphiasDef = data.frame("TSN" = unique(do.call("rbind", lapply(TSNsFromAphiasDef, "[[", 1))))
    TSNsFromAphiasDef$Aphia_ID = gsub("\\..*","",rownames(TSNsFromAphiasDef))
    names(TSNsFromAphiasDef)[names(TSNsFromAphiasDef) =="Aphia_ID"]<-"AphiaID"
    rownames(TSNsFromAphiasDef)<-NULL
    TSNsFromAphiasDef$TSN_src="worrms"
    TSNsFromAphiasDef$TSN_data="AphiaID"
    if (nrow(definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,])>0){
      definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,c("TSN","TSN_src","TSN_data")] = merge(definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,], TSNsFromAphiasDef, by.x="AphiaID", by.y="AphiaID")[,c("TSN.y","TSN_src.y","TSN_data.y")]
      definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,"TSN_Definitive"]<-FALSE
    }
  }
  #try to find remaining TSNs using sci and comm names
 
  definitive = checkAndCompareTSN(filterStrings, definitive, df=definitive, searchtype="scientific")
  missingTSN = definitive[!(definitive$TSN_Definitive %in% TRUE),]
  if (nrow(missingTSN)>0){
    definitive = checkAndCompareTSN(filterStrings, definitive, df=missingTSN, searchtype="common")
  }
  
  multi_Codes = multiCheck(multi_Codes)
  if (nrow(multi_Codes)>0){
    print("looking for the multimatch aphias")
    TSNsFromAphiasMulti = sapply_pb(multi_Codes$AphiaID,chk_worrmsTSN)
    TSNsFromAphiasMulti = data.frame("TSN" = unique(do.call("rbind", lapply(TSNsFromAphiasMulti, "[[", 1))))
    TSNsFromAphiasMulti$Aphia_ID = rownames(TSNsFromAphiasMulti)
    names(TSNsFromAphiasMulti)[names(TSNsFromAphiasMulti) =="Aphia_ID"]<-"AphiaID"
    rownames(TSNsFromAphiasMulti)<-NULL
    TSNsFromAphiasMulti$TSN_src="worrms"
    TSNsFromAphiasMulti$TSN_data="AphiaID"
    if (nrow(multi_Codes[multi_Codes$AphiaID %in% TSNsFromAphiasMulti$AphiaID,])>0){
      multi_Codes[multi_Codes$AphiaID %in% TSNsFromAphiasMulti$AphiaID,c("TSN","TSN_src","TSN_data")] = merge(multi_Codes[multi_Codes$AphiaID %in% TSNsFromAphiasMulti$AphiaID,], TSNsFromAphiasMulti, by.x="AphiaID", by.y="AphiaID")[,c("TSN.y","TSN_src.y","TSN_data.y")]
    }
  }
  
    multi_Codes= multiCheck(multi_Codes)
  if (nrow(multi_Codes)>0){
     cat("\nReview the contents of 'multi_Codes', as some matches were not definitive, and have alternates you should consider\n")
  }else{
    rm(multi_Codes)
  }
    
  mysterySpecCln = mysterySpec[mysterySpec$SPEC %in% (unique(c(definitive[definitive$TSN_Definitive %in% c(NA,FALSE),"SPEC"], definitive[definitive$AphiaID_Definitive %in% c(NA,FALSE),"SPEC"]))),"SPEC"]
  assign("mysterySpec", mysterySpecCln, envir = .GlobalEnv)
  names(definitive)[names(definitive)=="SPEC"]<-sci_col
  names(definitive)[names(definitive)=="COMM"]<-comm_col
  res = merge(spec_list_orig, definitive, by=c(sci_col,comm_col))
  return(res)
}  

