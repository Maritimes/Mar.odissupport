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
#' \item \\(.*?\\) -removes thing in brackets (e.g "(NS)")
#' \item \\b[a-zA-Z]{1,2}\\. - removes one or two letter blocks of text (potentially followed by a period) (e.g "SP.")
#' \item \\,\\s?(SMALL|LARGE) - removes instances like ",SMALL"
#' \item UNIDENTIFIED - removes the word "UNIDENTIFIED"
#' \item UNID\\. - removes the word "UNID."
#' \item EGGS - removes the word "EGSS"
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
  filterStrings <- c(filterStrings,"\\(.*?\\)", "(\\b[a-zA-Z]{1,2}\\.)","\\,\\s?(SMALL|LARGE)","-?UNIDENTIFIED", "\\,?\\s?EGGS-?","-?UNID\\.","PURSE\\s")
  # Initial data preparation ------------------------------------------------
  spec_list_orig <- spec_list  
  names(spec_list)[names(spec_list)==sci_col]<-"SPEC_ORIG"
  names(spec_list)[names(spec_list)==comm_col]<-"COMM_ORIG"
  spec_list$trackID <- seq.int(nrow(spec_list))
  spec_list$SPEC=trimws(gsub(paste(filterStrings,collapse="|"),  " ", spec_list$SPEC_ORIG))
  spec_list$COMM=trimws(gsub(paste(filterStrings,collapse="|"),  " ", spec_list$COMM_ORIG))

  definitive <- spec_list[,c("trackID","SPEC","COMM")]
  definitive<-definitive[!nchar(definitive$SPEC)==0,] 
  definitive$AphiaID_Definitive <- NA
  definitive$AphiaID <- NA
  definitive$AphiaID_src <- NA
  definitive$AphiaID_data <- NA
  definitive$TSN_Definitive <- NA
  definitive$TSN <- NA
  definitive$TSN_src <- NA
  definitive$TSN_data <- NA
 
  mysterySpec <-definitive
  resultFormat = data.frame(internal_name = character(), suggested_name = character(), search_type = character(),src= character(), ID= integer())
  
   multi_Codes <<- data.frame(suggested_name = character(), TSN= integer(), TSN_internal_name = character(), TSN_src= character(), TSN_data = character(),
                             AphiaID_internal_name = character(), AphiaID_src= character(), AphiaID_data = character(), AphiaID= integer() )
  # Basic helper functions --------------------------------------------------
  
  sapply_pb <- function(X, FUN, ...)
  {
    #this is a progress bar function I can wrap around stuff
    #think I stole it from https://ryouready.wordpress.com/2010/01/11/
    #progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/
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
  
  multiCheck<-function(multi_Codes, definitive, theCode){
if (nrow(multi_Codes)>0)browser()
    #update multi_codes to hold only those species for which we don't have a definitve value 
    theseMulti = multi_Codes[(multi_Codes$trackID %in% definitive[definitive[,theCode] %in% c(NA,FALSE),"trackID"]),]
    assign("multi_Codes", theseMulti, envir = .GlobalEnv)
    return(theseMulti)
  }
 
  # Functions for checking web services -------------------------------------
  chk_taxize <- function(info = NULL,searchtype='scientific', ask=FALSE, verbose=FALSE){
    info=strsplit(info,"\\$\\$\\$")
    searchterm=info[[1]][1]
    trackID=info[[1]][2]
    
    #use taxize functions to check stuff
    colname= ifelse(searchtype =="scientific","SPEC","COMM")
    this <- tryCatch(
      {
        taxize::get_wormsid(query=searchterm, searchtype=searchtype, ask=ask, verbose=verbose, accepted = TRUE,rows = 1, messages =FALSE)
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
        #mult$trackID = trackID
        mult$AphiaID = ID
        mult$AphiaID_src = "Taxize"
        mult$AphiaID_data = searchtype
        mult$AphiaID_internal_name = searchterm
        mult$suggested_name = NA #can't get names back
        mult$TSN_internal_name = NA
        mult$TSN_src = NA
        mult$TSN_data = NA
        mult$TSN = NA
        mult$SPEC<-NULL
        mult$COMM<-NULL
        names(mult)[names(mult) == 'SPEC_ORIG'] <- 'SPEC'
        names(mult)[names(mult) == 'COMM_ORIG'] <- 'COMM'
        multi_Codes <- rbind(.GlobalEnv$multi_Codes, mult)
        assign("multi_Codes", multi_Codes, envir = .GlobalEnv)
        
        #can't capture all other matches - only can get first
      } else if (attr(this,"match")=="found" & attr(this,"multiple_matches") == FALSE) {
        cmt = "multi match - kept only accepted"
      }else{
        cmt = "multi match - no ids returned"
      }
    }
    this = as.data.frame(cbind("SEARCHTERM" = searchterm, "ID_taxize"=ID,"MATCH_taxize"= cmt, trackID= trackID))
    if (NROW(this[grep("NA due to ask=FALSE", this$MATCH_WORMS),])>0)  this[grep("NA due to ask=FALSE", this$MATCH_WORMS),]$MATCH_WORMS<-"multi match"
    return(this)
  }
  
  chk_worrms <- function(info = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    info=strsplit(info,"\\$\\$\\$")
    searchterm=info[[1]][1]
    trackID=info[[1]][2]
    colname= ifelse(searchtype =="scientific","SPEC","COMM")
    resName = "scientificname"
    if (searchtype == 'scientific'){
      thisi <- tryCatch(
        {
          worrms::wm_records_name(searchterm, fuzzy=T, marine_only = T)
        },
        error=function(cond){
        }
      )
    }else{
      thisi <- tryCatch(
        {
          worrms::wm_records_common(searchterm, fuzzy=T,marine_only = T)
        },
        error=function(cond){
        }
      )
    }
    crapRec = data.frame("SEARCHTERM" = searchterm, "ID_worrms" = NA, "MATCH_worrms" = "not found","trackID"= trackID)
    if (is.null(thisi)) {
      thisi <- crapRec
    } else if (nrow(thisi)==0) {
      thisi <- crapRec
    }else if (nrow(thisi)==1) {
      thisi = cbind("SEARCHTERM"= searchterm,"ID_worrms" = thisi$valid_AphiaID, "MATCH_worrms" = "single match","trackID"= trackID)
    }else{
      if (nrow(thisi[thisi$status=="accepted",])==1){
        #maybe need merge thisi with speclist?
        thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi[thisi$status=="accepted","valid_AphiaID"], "MATCH_worrms" = "multi match - kept only accepted","trackID"= trackID)
        colnames(thisi)[colnames(thisi)=="AphiaID"] <- "ID_worrms"
      } else {
        mult = thisi[,c("scientificname","valid_AphiaID")]
        colnames(mult)<-c("suggested_name","AphiaID")
        mult$AphiaID_src = "worrms"
        mult$AphiaID_internal_name = searchterm
        mult$AphiaID_data = searchtype
        mult$TSN_internal_name = NA
        mult$TSN_src = NA
        mult$TSN_data = NA
        mult$trackID = trackID
        # mult$TSN = NA
        mult = merge(spec_list, mult, by.x="trackID", by.y="trackID")
        mult$SPEC<-NULL
        mult$COMM<-NULL
        names(mult)[names(mult) == 'SPEC_ORIG'] <- 'SPEC'
        names(mult)[names(mult) == 'COMM_ORIG'] <- 'COMM'
        multi_Codes <- rbind(.GlobalEnv$multi_Codes, mult)
        assign("multi_Codes", multi_Codes, envir = .GlobalEnv)
        if (nrow(thisi[thisi$status=="accepted",])>1) {
          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - kept 1st accepted","trackID"= trackID)
        }else {
          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - none accepted; retained first","trackID"= trackID)
        }
      }
    }
    return(thisi)
  }
  
  chk_ritis <- function(info = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    info=strsplit(info,"\\$\\$\\$")
    searchterm=info[[1]][1]
    trackID=info[[1]][2]
    
    this = resultFormat
    if (searchtype == 'scientific'){
      resName = "combinedName"
      colname = "SPEC"
      thisi =  as.data.frame(ritis::search_scientific(searchterm, wt = "json", raw = FALSE))
    }else{
      resName = "commonName"
      colname = "COMM"
      if (searchterm=="WHALES")browser()
      thisi = as.data.frame(ritis::search_common(searchterm, wt = "json", raw = FALSE))
    }
    if (nrow(thisi)==0) {
      thisi=data.frame("SEARCHTERM" = searchterm, "ID_ritis" = NA, "MATCH_ritis" = "not found","trackID"= trackID)
    }else if (nrow(thisi)==1) {
      thisi = cbind("SEARCHTERM"= searchterm,"ID_ritis" = thisi$tsn, "MATCH_ritis" = "single match","trackID"= trackID)
    }else{
      thisi[,resName]<-toupper(thisi[,resName])
      mult = thisi[,c(resName,"tsn")]
      colnames(mult)<-c("suggested_name","TSN")
      mult$TSN_internal_name = searchterm
      #mult$joiner = searchterm
      mult$TSN_src = "ritis"
      mult$TSN_data = searchtype
      mult$AphiaID_internal_name = NA
      mult$AphiaID_src = NA
      mult$AphiaID_data = NA
      mult$AphiaID = NA
      mult$trackID = trackID
      mult = merge(spec_list, mult, by.x="trackID", by.y="trackID")
      mult$SPEC<-NULL
      mult$COMM<-NULL
      names(mult)[names(mult) == 'SPEC_ORIG'] <- 'SPEC'
      names(mult)[names(mult) == 'COMM_ORIG'] <- 'COMM'
      mult$TSN<-mult$TSN.x
      mult$TSN.x<-NULL
      mult$TSN.y<-NULL
      #multi_Codes <<- rbind(multi_Codes, mult)
      if (length(thisi[thisi[,resName] == searchterm,"tsn"])==1){
        thisi = cbind("SEARCHTERM"=searchterm, "ID_ritis" = thisi[thisi[,resName] == searchterm,"tsn"], "MATCH_ritis" ="multi match - retained identical spelling","trackID"= trackID)
      }else{
        thisi = cbind("SEARCHTERM"=searchterm, "ID_ritis" =thisi$tsn[1], "MATCH_ritis" ="multi match - no identical spelling; retained first","trackID"= trackID)
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
    colname= ifelse(searchtype =="scientific","SPEC","COMM")
    searchtypeVal= ifelse(colname =="SPEC","scientific name","common name")

    if(!any(grep("trackID",colnames(curDefinitive)))){
      curDefinitive$searchterm<-NA
      curDefinitive$trackID<-NA
      for (i in 1:nrow(curDefinitive)){
      curDefinitive[i,"searchterm"]=strsplit(curDefinitive[,colname],"\\$\\$\\$")[[i]][1]
      curDefinitive[i,"trackID"]=strsplit(curDefinitive[,colname],"\\$\\$\\$")[[i]][2]
      }
    }
    if(any(grep("taxize",colnames(results))) |(any(grep("worrms",colnames(results))) )){
      theID = "AphiaID"
      theIDComment = "AphiaID_src"
      theIDUsed = "AphiaID_data"
      isDefFlag = "AphiaID_Definitive"
      otherTheID = "TSN"
      otherTheIDComment = "TSN_src"
      otherTheIDUsed = "TSN_data"
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
      isDefFlag = "TSN_Definitive"
      otherTheID = "AphiaID"
      otherTheIDComment = "AphiaID_src"
      otherTheIDUsed = "AphiaID_data"
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
    results <-results[,c("SEARCHTERM","trackID", theID,theIDComment,theIDUsed, isDefFlag)]
    names(results)[names(results) == 'SEARCHTERM'] <- colname
    
    #multimatches = 
    #did something get deleted here?!
 
    #just adding new stuff - we don't have anything yet

    newFind = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(NA)),][c("trackID",theID)], 
                    results[results[,isDefFlag] %in% c(TRUE,FALSE),], by="trackID")
    newFind[,'thisID'] <- newFind[paste0(theID,".y")]
     
    #if the new codes are definitive, we will replace the existing (not-definitive) ones entirely 

    replacers = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(FALSE)),][c("trackID",theID)], 
                      results[results[,isDefFlag] %in% TRUE,], by="trackID")
    replacers[,'thisID'] <- replacers[paste0(theID,".y")]
          
    #if we have already have a tentative code, and our new result set gives the same code,
    #we will call it definitive   

    doubleblind = merge(curDefinitive[(curDefinitive[,isDefFlag] %in% c(FALSE)),][c("trackID",theID)], 
                        results[results[,isDefFlag] %in% FALSE,], by="trackID")
    doubleblind[,'thisID'] <- doubleblind[paste0(theID,".y")]    
     
    
    newCurDef<-curDefinitive
    repCurDef<-curDefinitive
    dbCurDef<-curDefinitive
    

    if (nrow(newFind)>0){
      newCurDef = merge(newCurDef[,!(names(newCurDef) %in% c(theID,theIDComment,theIDUsed, isDefFlag))], 
                        newFind[,(names(newFind) %in% c("thisID","trackID",theIDComment,theIDUsed, isDefFlag))], by="trackID")
      names(newCurDef)[names(newCurDef) == 'thisID'] <- theID
    }else{
      newCurDef <- newCurDef[0,]
    }
    if (nrow(replacers)>0){
      repCurDef = merge(repCurDef[,!(names(repCurDef) %in% c(theID,theIDComment,theIDUsed, isDefFlag))], 
                        replacers[,(names(replacers) %in% c("thisID","trackID",theIDComment,theIDUsed, isDefFlag))],  by="trackID")
      names(repCurDef)[names(repCurDef) == 'thisID'] <- theID
    }else{
      repCurDef <- repCurDef[0,]
    }
    
    if (nrow(doubleblind)>0){
      dbCurDef = merge(dbCurDef[,!(names(dbCurDef) %in% c(theID,theIDUsed, isDefFlag))], 
                       doubleblind[,(names(doubleblind) %in% c("thisID","trackID",theIDComment,theIDUsed, isDefFlag))],by="trackID")
      commentx = paste0(theIDComment,".x")
      commenty = paste0(theIDComment,".y") 
      
      dbCurDef[,theIDComment]=paste0(dbCurDef[,commentx], " & ", dbCurDef[,commenty])
      dbCurDef[,commentx]<-dbCurDef[,commenty]<-NULL
      names(dbCurDef)[names(dbCurDef) == 'thisID'] <- theID
      dbCurDef[,isDefFlag] <-TRUE
    } else{
      dbCurDef <- dbCurDef[0,]
    }

    newStuff<-rbind(newCurDef, repCurDef, dbCurDef)
    curDefinitive=rbind(curDefinitive[!(curDefinitive$trackID %in% newStuff$trackID),],newStuff)
    return(curDefinitive)
  }   
  checkAndCompareTSN<-function(definitive=definitive, df=NULL, searchtype = NULL){
    colname= ifelse(searchtype =="scientific","SPEC","COMM")
    print(paste("...Checking ritis using ",searchtype," names"))
    no_TSN = definitive[definitive$TSN_Definitive %in% c(NA, FALSE),]
    if (nrow(no_TSN)>0) {
      info=within(no_TSN, x <- paste(get(colname),trackID,sep='$$$'))
      no_TSN = as.data.frame(t(sapply_pb(info$x,chk_ritis, searchtype=searchtype)))
      names(no_TSN) = c("SEARCHTERM","ID_ritis","MATCH_ritis","trackID")
      definitive = assignDefinitive(definitive,no_TSN, searchtype)
    }
    return(definitive)
  }
  
  checkAndCompareAphiaID<-function(definitive, df= NULL, searchtype = NULL){
    colname= ifelse(searchtype =="scientific","SPEC","COMM")
    no_aphiaID = df[!(df$AphiaID_Definitive %in% TRUE),]
    if (nrow(no_aphiaID)>0) {
      print(paste0("...against taxize using ",searchtype," names"))
      info=within(no_aphiaID, x <- paste(get(colname),trackID,sep='$$$'))
      no_aphiaID = as.data.frame(t(sapply_pb(info$x,chk_taxize, searchtype=searchtype)))
    }
    if (nrow(no_aphiaID)>0) {
      no_aphiaID = data.frame("SEARCHTERM" = toupper(unlist(no_aphiaID$SEARCHTERM)),
                              "ID_taxize" = unlist(no_aphiaID$ID_taxize),
                              "MATCH_taxize" = unlist(no_aphiaID$MATCH_taxize),
                              "trackID" = unlist(no_aphiaID$trackID))
      definitive = assignDefinitive(definitive,no_aphiaID, searchtype)
      not_definitive = definitive[!(definitive$AphiaID_Definitive %in% TRUE),]
      if (nrow(not_definitive)>0){
        print(paste0("...against worrms using ",searchtype," names"))
        info=within(not_definitive, x <- paste(get(colname),trackID,sep='$$$'))
        not_definitive = as.data.frame(t(sapply_pb(info$x,chk_worrms, searchtype=searchtype)))
        names(not_definitive)<-c("SEARCHTERM", "ID_worrms", "MATCH_worrms", "trackID")
        not_definitive = data.frame("SEARCHTERM" = toupper(unlist(not_definitive$SEARCHTERM)),
                                    "ID_worrms" = unlist(not_definitive$ID_worrms),
                                    "MATCH_worrms" = unlist(not_definitive$MATCH_worrms),
                                    "trackID" = unlist(not_definitive$trackID))
        definitive = assignDefinitive(definitive,not_definitive, searchtype)
      }
    }
    return(definitive)
  }  
  
  # Start of the actual processing ------------------------------------------

  print("Looking for AphiaIDs")
  definitive = checkAndCompareAphiaID(definitive, df = definitive,searchtype="scientific")
  missing = definitive[!(definitive$AphiaID_Definitive %in% TRUE),]

  multi_Codes = multiCheck(multi_Codes, definitive, "AphiaID")
  if (nrow(missing)>0){
    definitive = checkAndCompareAphiaID(definitive, df = missing,searchtype="common")
    multi_Codes = multiCheck(multi_Codes, definitive, "AphiaID")
  }
definitive_postAphia<<-definitive  
  # #have lots of aphias - try to use them to get TSNs
  print("Looking for TSNs")
  print("...Checking worrms using AphiaIDs")
  if (!is.null(definitive)){
    TSNsFromAphiasDef = sapply_pb(definitive$AphiaID,chk_worrmsTSN)
    TSNsFromAphiasDef = data.frame("TSN" = unique(do.call("rbind", lapply(TSNsFromAphiasDef, "[[", 1))))

    TSNsFromAphiasDef$Aphia_ID = gsub("\\..*","",rownames(TSNsFromAphiasDef))
    names(TSNsFromAphiasDef)[names(TSNsFromAphiasDef) =="Aphia_ID"]<-"AphiaID"
    rownames(TSNsFromAphiasDef)<-NULL
    if (nrow(TSNsFromAphiasDef)==0){

    }else{
      TSNsFromAphiasDef$TSN_src="worrms"
      TSNsFromAphiasDef$TSN_data="AphiaID"
      if (nrow(definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,])>0){
        definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,c("TSN","TSN_src","TSN_data")] = merge(definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,], TSNsFromAphiasDef, by.x="AphiaID", by.y="AphiaID")[,c("TSN.y","TSN_src.y","TSN_data.y")]
        definitive[definitive$AphiaID %in% TSNsFromAphiasDef$AphiaID,"TSN_Definitive"]<-FALSE
      }
    }
  }
  #try to find remaining TSNs using sci and comm names
  definitive = as.data.frame(sapply(checkAndCompareTSN(definitive, df=definitive, searchtype="scientific"),unlist))
  missingTSN = definitive[!(definitive$TSN_Definitive %in% TRUE),]
  if (nrow(missingTSN)>0){
    info=within(missingTSN, x <- paste(COMM,trackID,sep='$$$'))
    definitive = as.data.frame(sapply(checkAndCompareTSN(definitive, df=info, searchtype="common"),unlist))
  }
  multi_Codes = multiCheck(multi_Codes, definitive, "TSN")
  
  definitive_postTSN<<-definitive  
  if (nrow(multi_Codes)>0){
    
    multi_Codes_postAll<-multi_Codes
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
    multi_Codes= multiCheck(multi_Codes, definitive,"TSN")
  }

  if (nrow(multi_Codes)>0){
     cat("\nReview the contents of 'multi_Codes', as some matches were not definitive, and have alternates you should consider\n")
  }else{
    rm(multi_Codes)
  }
  mysterySpecCln = mysterySpec[mysterySpec$SPEC %in% (unique(c(definitive[definitive$TSN_Definitive %in% c(NA,FALSE),"SPEC"], definitive[definitive$AphiaID_Definitive %in% c(NA,FALSE),"SPEC"]))),"SPEC"]
  assign("mysterySpec", mysterySpecCln, envir = .GlobalEnv)
  
  names(definitive)[names(definitive)=="SPEC"]<-sci_col
  names(definitive)[names(definitive)=="COMM"]<-comm_col
  res = merge(spec_list[,!(colnames(spec_list) %in% c("TSN", "COMM", "SPEC"))], definitive, by="trackID")
  res$SPEC<-NULL
  res$COMM<-NULL
  names(res)[names(res) == 'SPEC_ORIG'] <- 'SPEC'
  names(res)[names(res) == 'COMM_ORIG'] <- 'COMM'
  
  return(res)
}  

