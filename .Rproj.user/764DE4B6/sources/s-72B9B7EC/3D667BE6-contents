#' @title getTaxaIDs
#' @description This function attempts to determine the aphiaIDs (WoRMS) and
#' TSNs (ITIS) using first the scientific name, and then the common name. If an
#' aphiaid is found, but no TSN, it will try to find the TSN using the aphiaid.
#' If a TSN is found, but no aphiaid, it will try to find the aphiaid using the
#' TSN.
#' Additionally, in addition to all of the original fields, the returned data
#' will also include the WoRMS aphiaid, the ITIS TSN, and will indicate the
#' method used to find each.  Also for each, a field will indicate whether one,
#' many, or no matches were found.  Lastly, if the species has a different
#' "accepted name" than the one provided, that will also be returned.
#' @param spec_list the dataframe containing information to be decoded to TSN and
#' aphiaIDs
#' @param sci_col the name of the column of the dataframe containing the
#' scientific names
#' @param comm_col the name of the column of the dataframe containing the
#' common names
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
getTaxaIDs <- function(spec_list=NULL, sci_col=NULL, comm_col=NULL, debug = FALSE){
  #' remove potential confusion from CASE sensitivity
  if (nrow(spec_list)==1){
    spec_list = as.data.frame(t(as.data.frame(sapply(spec_list, toupper))))
  }else{
  spec_list = as.data.frame(sapply(spec_list, toupper))
}
  resultFormat = data.frame(internal_name = character(), suggested_name = character(), search_type = character(),src= character(), ID= integer())

  #' "definitive_Aphia" will be our final results where we are confident about our 
  #' results.  These will include where the series of checks coalesce on a 
  #' single "accepted" value if one check (e.g. scientific name) finds multiple 
  #' (or no) ids, another check will be performed (eg common name).  If the 
  #' follow up check finds a single match, the result will be added to the 
  #' definitive_Aphia list
  #' 
  #' "multiMatches" will hold IDs for which we never found a definitive_Aphia value.  
  #' For example, if multiple IDS were found in all checks, and no one of them 
  #' was deemed "accepted", all of these options will be added to multiMatches.  
  #' Similarly, if all of the checks result in a group of "accepted matches" 
  #' that are different from one another, the results will also be stored.  The
  #' expectation is that a human should go through the list of options and 
  #' pick the best.
  #' 
#'_src is where we got the ID (e.g. worrms)
#'_data is what we used to get the id (e.g. scientific/common name, or another ID)
  definitive_Aphia <<- data.frame(SPEC = character(), AphiaID = integer(), AphiaID_src = character(), AphiaID_data = character())
  definitive_TSN <<- data.frame(SPEC = character(), TSN = integer() ,TSN_src = character(), TSN_data = character())
  
  
  
  #multimatches will hold alternative id values for species where we found multiple matches
  multi_Codes <<- data.frame(suggested_name = character(), TSN= integer(), TSN_internal_name = character(), TSN_src= character(), TSN_data = character(),
                             AphiaID_internal_name = character(), AphiaID_src= character(), AphiaID_data = character(), AphiaID= integer() )

  mysterySpec <<-spec_list
  
  if (is.null(comm_col)){
      df_imp<-data.frame("tmpzzz" = spec_list[,sci_col])
      names(df_imp)[names(df_imp)=="tmpzzz"]<-sci_col
  } else {
    df_imp<-data.frame(spec_list[,c(sci_col, comm_col)])
  }

  definitive <<- df_imp
  definitive$AphiaID <- NA
  definitive$AphiaID_src <- NA
  definitive$AphiaID_data <- NA
  definitive$TSN <- NA
  definitive$TSN_src <- NA
  definitive$TSN_data <- NA
  assign("definitive", definitive, envir = .GlobalEnv)
  df_res<-spec_list

  sapply_pb <- function(X, FUN, ...)
  {
    #this is a progress bar function I can wrap around stuff
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
  
  # Functions for checking web services -------------------------------------
  chk_taxize <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    
    if (debug) cat("\nin chk_taxize")
    #use taxize functions to check stuff
    if (grepl("\\(NS\\)", searchterm)==TRUE){
      tmpTerm = trimws(gsub("\\(NS\\)",  "", searchterm))
    }else{
      tmpTerm = searchterm
    }
    this <- tryCatch(
      {
        taxize::get_wormsid(query=tmpTerm, searchtype=searchtype, ask=ask, verbose=verbose, accepted = TRUE,rows = 1)
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
        #can't capture all other matches - only can get first
      } else if (attr(this,"match")=="found" & attr(this,"multiple_matches") == FALSE) {
        cmt = "multi match - kept only accepted"
      }else{
        cmt = "multi match - no ids returned"
      }
    }
    this = as.data.frame(cbind("SEARCHTERM" = searchterm, "ID_taxize"=ID,"MATCH_taxize"= cmt))
    if (NROW(this[grep("NA due to ask=FALSE", this$MATCH_WORMS),])>0)  this[grep("NA due to ask=FALSE", this$MATCH_WORMS),]$MATCH_WORMS<-"multi match"
    
    if (debug) cat("\ndone chk_taxize")
    return(this)
  }
  chk_worrms <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #browser()
    #use worrms functions to check stuff
    if (debug) cat("\nin chk_worrms")
    #this = resultFormat
    if (searchtype == 'scientific'){
      colname = "SPEC"
    }else{
      colname = "COMM"
    }
    # 
    # for (i in 1:NROW(searchterm)){
    #   if (grepl("\\(NS\\)", searchterm[i])==TRUE){
    #     tmpTerm = trimws(gsub("\\(NS\\)",  "", searchterm[i]))
    #   }else{
    #     tmpTerm = searchterm[i]
    #   }
       if (grepl("\\(NS\\)", searchterm)==TRUE){
         tmpTerm = trimws(gsub("\\(NS\\)",  "", searchterm))
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
            #message(cond)
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
          #thisi = merge(thisi, spec_list, by.x="SEARCHTERM", by.y=colname )
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
                mult = merge(spec_list, mult, by.x=colname, by.y="joiner")
                mult$joiner = NULL
                mult[,colname] = mult$TSN_internal_name
                multi_Codes <<- rbind(multi_Codes, mult)
                
                        if (nrow(thisi[thisi$status=="accepted",])>1) {
                          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - kept 1st accepted")
                        }else {
                          thisi = cbind("SEARCHTERM" = searchterm,"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - none accepted; retained first")
                        }
              }
        }
    #}
      if (debug) cat("\ndone chk_worrms")
    return(thisi)
  }
  
  chk_ritis <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use ritis functions to check stuff
    if (debug) cat("\nin chk_ritis")
    this = resultFormat
    # for (i in 1:NROW(searchterm)){
      
      if (grepl("\\(NS\\)", searchterm)==TRUE){
        tmpTerm = trimws(gsub("\\(NS\\)",  "", searchterm))
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
      #make sure we're looking at a common case
      
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
        mult$joiner = NULL
        multi_Codes <<- rbind(multi_Codes, mult)
        if (length(thisi[thisi[,resName] == searchterm,"tsn"])==1){
          thisi = cbind(searchterm, thisi[thisi[,resName] == searchterm,"tsn"], "multi match - retained identical spelling")
        }else{
          thisi = cbind(searchterm, thisi$tsn[1], "multi match - no identical spelling; retained first")
        }
      }
    #   this = rbind(this, thisi)
    # }
    #browser()
    #this = this[rowSums(is.na(this)) != ncol(this),]
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
    if (debug) cat("\ndone chk_ritis")
    return(thisi)
  }
  
  
  assignDefinitive<-function(results = NULL, searchtype = NULL){
    if (debug) cat("\nin assignDefinitive")
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
      defTableName = "definitive_Aphia"
      
      otherTheID = "TSN"
      otherTheIDComment = "TSN_src"
      otherTheIDUsed = "TSN_data"
      otherDefTableName = "definitive_TSN"
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
        #if (colname =="COMM") browser()
      theID = "TSN"
      theIDComment = "TSN_src"
      theIDUsed = "TSN_data"
      defTableName = "definitive_TSN"
      otherTheID = "AphiaID"
      otherTheIDComment = "AphiaID_src"
      otherTheIDUsed = "AphiaID_data"
      otherDefTableName = "definitive_Aphia"
        
      results_ID = "ID_ritis"
      results_MATCH = "MATCH_ritis"
      results_type = "ritis"
    }
    
    excludes <- c(theID,theIDComment,theIDUsed)   
    keeperStrings = c("single match","only accepted")
    # results = merge(df1, df2, by="SEARCHTERM")
    results[theID] = NA
    results[theIDComment] = NA
    results[grep(paste(keeperStrings,collapse="|"),results[,results_MATCH]),theID]<-results[grep(paste(keeperStrings,collapse="|"),results[,results_MATCH]),results_ID]
    
    #no ID = unmatched
    if(nrow(results[is.na(results[,results_ID]),])>0) results[is.na(results[,results_ID]),theIDComment]<-"Unmatched"
    #only 1 id present, and it was in the keeper strings (i.e. a confident match)

    if(nrow(results[(!is.na(results[,results_ID]) & !is.na(results[,theID])),])>0)
      results[(!is.na(results[,results_ID]) & !is.na(results[,theID])),theIDComment] <-paste(results_type)
    #all that remain should be those which had 1 or more matches, but none within the keeper strings
    
    if(nrow(results[is.na(results[,theIDComment]),])>0)
      results[is.na(results[,theIDComment]),theIDComment]<-"Ambiguous Match"
    
    results[,theIDUsed] =NA

     #check the use of theID, below
    if (nrow(results[!is.na(results[theID]),])>0) results[!is.na(results[theID]),theIDUsed]<-searchtypeVal
    #results_MATCH
    results <-results[,c("SEARCHTERM",theID,theIDComment,theIDUsed)]
    names(results)[names(results) == 'SEARCHTERM'] <- colname
    #'get the records we have confident IDs for and add them to the appropriate 
    #'definitive table (which is a global object)
    thesedefinitive = results[!(results[,theIDComment] %in% c("Ambiguous Match", "Unmatched")),]
    #and save the others
    theseunhandled = results[!(results[,colname] %in% thesedefinitive[,colname]),]
    definitive =.GlobalEnv$definitive

    #newdefinitive = definitive
    cln_definitive = definitive[is.na(definitive[,theID]),]
    keep_definitive = definitive[!is.na(definitive[,theID]),]
  
    if (nrow(thesedefinitive)>0){
  
      #only get NA ids so we don't overwrite values
      #merge the remaining definitive values with the newly found ones(skipping the columns we will be populating)
      new_definitive = merge(cln_definitive[!names(cln_definitive) %in% excludes], thesedefinitive, by=colname)
      new_definitive = new_definitive[,c("SPEC",excludes)]
      
      definitive[definitive$SPEC %in% new_definitive$SPEC, excludes]<-new_definitive[excludes]
      #write the new values to the original data frame
      # newdefinitive[newdefinitive[,colname] %in% cln_definitive[,colname],excludes]<-cln_definitive[,excludes]
      # newdefinitive <- newdefinitive[!is.na(newdefinitive[,theID]),]
       
      assign("definitive", definitive, envir = .GlobalEnv )
    }  
    if (debug) cat("\ndone assignDefinitive")
    return(theseunhandled)
  }   

  checkAndCompare<-function(df= NULL, searchtype = NULL){
    
    if (debug) cat("\nin checkAndCompare")
    if (searchtype =="scientific"){
      colname = "SPEC"
    }else{
      colname = "COMM"
    }
df=merge(df,spec_list)
    # check for matches based on sci name, then common name (if provided)
    cat(paste("Checking using ",searchtype," names\n"))
    cat("\tagainst Taxize (AphiaID)....\n")
    no_ID = as.data.frame(t(sapply_pb(df[,colname],chk_taxize, searchtype=searchtype)))
    no_ID = data.frame("SEARCHTERM" = toupper(unlist(no_ID$SEARCHTERM)),
                                  "ID_taxize" = unlist(no_ID$ID_taxize),
                                  "MATCH_taxize" = unlist(no_ID$MATCH_taxize))
    if (nrow(no_ID)>0) no_ID = assignDefinitive(no_ID, searchtype)
    #only look for those we're not sure of
    if (nrow(no_ID)>0) {
      cat("\tagainst worrms (AphiaID)\n")
      no_ID = as.data.frame(t(sapply_pb(no_ID[,colname],chk_worrms, searchtype=searchtype)))

      no_ID = data.frame("SEARCHTERM" = toupper(unlist(no_ID$SEARCHTERM)),
                                    "ID_worrms" = unlist(no_ID$ID_worrms),
                                    "MATCH_worrms" = unlist(no_ID$MATCH_worrms))
      if (nrow(no_ID)>0) no_ID = assignDefinitive(no_ID, searchtype)
    }
    if (nrow(no_ID)>0) {
      cat("\tagainst ritis (TSN) \n")
      

      no_ID = as.data.frame(t(sapply_pb(no_ID[,colname],chk_ritis, searchtype=searchtype)))
      
      
      if (searchtype == "scientific"){
      no_ID = data.frame("SEARCHTERM" = toupper(unlist(no_ID$SEARCHTERM)),
                         "ID_ritis" = unlist(no_ID$ID_ritis),
                         "MATCH_ritis" = unlist(no_ID$MATCH_ritis))
      }else{
        #not sure why the column names are different for sci vs comm
      no_ID = data.frame("SEARCHTERM" = toupper(unlist(no_ID$V1)),
                         "ID_ritis" = unlist(no_ID$V2),
                         "MATCH_ritis" = unlist(no_ID$V3))
      } 
      if (nrow(no_ID)>0) no_ID = assignDefinitive(no_ID, searchtype)
    }
    
    if (debug) cat("\ndone checkAndCompare")
    if (nrow(no_ID)==0){
      cat("\nAll species identified\n")
      return(NULL)
    } else {
      return(no_ID)
    }
  }

  multiCheck<-function(){
    allMulti = merge(spec_list, multi_Codes)
    allMulti = allMulti[!(allMulti$SPEC %in% definitive_Aphia$SPEC) | !(allMulti$SPEC %in% definitive_TSN$SPEC) ,]
    assign("multi_Codes",allMulti, envir = .GlobalEnv )
  }
# Start of actual check (all functions above) -----------------------------
    #'functions are going to populate some global variables - including:
    #'definitive_Aphia (aphias for which there is no confusion)
    #'definitive_TSN (TSNs for which there is no confusion)
    #'multi_Aphia (aphias for which there is confusion)
    #'
    #'multi_TSN (aphias for which there is confusion)
    
  #Check for aphias by scientific name

  missing = checkAndCompare(df = df_imp,searchtype="scientific")
  multiCheck()
  if (!is.null(missing)){
    cat("\nChecking ambiguously and unmatched records using Common names\n")
    missing = checkAndCompare(df = missing,searchtype="common")
    multiCheck()
  }
  #' Done our initial check for Aphias, now use the ones we found to find TSNs
  if (!is.null(missing)){  
    
    browser()
    cat("\nFinding TSNs using the 'definitive' AphiaIDs\n")
    TSNsFromAphias = data.frame(t(data.frame(sapply_pb(.GlobalEnv$definitive$AphiaID,chk_worrmsTSN)))[,1])
    TSNsFromAphias$AphiaID<-gsub("\\..*","",rownames(TSNsFromAphias)) 
    rownames(TSNsFromAphias)<-NULL
    colnames(TSNsFromAphias)<-c("TSN", "AphiaID")
    TSNsFromAphias$AphiaID<-gsub("X","",TSNsFromAphias$AphiaID) 
    TSNsFromAphias = unique(TSNsFromAphias)
    #TSNsFromAphias= data.frame(AphiaID = gsub("X","",rownames(TSNsFromAphias)),TSN = do.call(pmax, c(TSNsFromAphias, list(na.rm=TRUE))))
    theseDefinitive_TSN = data.frame(SPEC = definitive_Aphia$SPEC, 
                                TSN_src = "worrms", 
                                TSN_data = "AphiaID (Definitive)", 
                                AphiaID = definitive_Aphia$AphiaID, 
                                row.names = NULL)
  #some definitiveAphias are duplicated with only a different aphiaid_src
    theseDefinitive_TSN = merge(unique(theseDefinitive_TSN), TSNsFromAphias)
    theseDefinitive_TSN$AphiaID<-NULL
    theseDefinitive_TSN = theseDefinitive_TSN[!(theseDefinitive_TSN$SPEC %in% definitive_TSN$SPEC),]
    if (nrow(theseDefinitive_TSN)>0)definitive_TSN = rbind(definitive_TSN,theseDefinitive_TSN)
      assign("definitive_TSN",definitive_TSN, envir = .GlobalEnv )
      multiCheck()
    }

#'   
#'   cat("\nFinding TSNs using the AphiaIDs of multimatches\n") 
#'   #' now lets see if we can find TSNs for the multimatches
#' 
#'   multi_TSN = sapply_pb(multi_Codes$AphiaID,chk_worrmsTSN) 
#'   #data.frame(mm=names(TSNsFromAphias))
#'   multi_TSN = data.frame(TSN = unlist(lapply(multi_TSN, `[`,1)))
#'   multi_TSN$AphiaID<-gsub(".*\\.","",rownames(multi_TSN)) 
#'   rownames(multi_TSN)<-NULL
#'   multi_TSN$TSN_src = "worrms"
#'   multi_TSN$TSN_data = "AphiaID (multi)"
#'   multi_Codes = merge(multi_Codes, multi_TSN, all.x = T)
#'   
#'   #' get the records for which we got TSNs, and if there's only one, keep the 
#'   #' associated AphiaID
#'   multi_solved = count(multi_Codes[!is.na(multi_Codes$TSN),], "AphiaID_internal_name")
#'   multi_solved = multi_solved[multi_solved$freq ==1,"AphiaID_internal_name"]
#'   multi_solved = multi_Codes[multi_Codes$AphiaID_internal_name %in% multi_solved & !is.na(multi_Codes$TSN),]
#'   #if we're adding to the aphiaId, it's because the TSN was found, if we're 
#'   #adding to the TSN it's because we found it from the AphiaID
#'   multi_solved$AphiaID_data = "TSN (MultiMatch)"
#'   multi_solved$TSN_data = "AphiaID (MultiMatch)"
#'   
#'   names(multi_solved)[names(multi_solved) == 'AphiaID_internal_name'] <- 'SPEC'
#'   definitive_Aphia = rbind(definitive_Aphia, multi_solved[,c("SPEC", "AphiaID","AphiaID_src","AphiaID_data")])
#'   assign("definitive_Aphia",definitive_Aphia, envir = .GlobalEnv)
#'   
#'   
#'   names(multi_solved)[names(multi_solved)=="AphiaID_src"]<-"TSN_src"
#'   definitive_TSN = rbind(definitive_TSN, multi_solved[,c("SPEC", "TSN","TSN_src","TSN_data")])
#'   assign("definitive_TSN",definitive_TSN, envir = .GlobalEnv)
#' 
#'   
#' # Now look for TSNs for unmatched spec using Sci and Comm names -----------
#' 
#'   remaining = df_imp[!(df_imp$SPEC %in% definitive_TSN$SPEC | df_imp$SPEC %in% definitive_Aphia$SPEC), ]
#'   cat("\nChecking scientific names of remaining unmatched species ritis (TSN) \n")
#'       spec_list_ritis_sci = as.data.frame(t(sapply_pb(remaining[,sci_col],chk_ritis)))
#'       spec_list_ritis_sci = data.frame("SEARCHTERM" = toupper(unlist(spec_list_ritis_sci$V1)),
#'                                        "ID_ritis" = unlist(spec_list_ritis_sci$V2),
#'                                        "MATCH_ritis" = unlist(spec_list_ritis_sci$V3))
#'       remaining1=remaining[!(remaining$SPEC %in% spec_list_ritis_sci$SEARCHTERM),]
#'   cat("\nChecking common names of remaining unmatched species ritis (TSN) \n")
#'        spec_list_ritis_comm = as.data.frame(t(sapply_pb(remaining1[,comm_col],chk_ritis, searchtype="common")))
#'        spec_list_ritis_comm = data.frame("SEARCHTERM" = toupper(unlist(spec_list_ritis_comm$V1)),
#'                                         "ID_ritis" = unlist(spec_list_ritis_comm$V2),
#'                                         "MATCH_ritis" = unlist(spec_list_ritis_comm$V3))
#'        

       #spec_list_ritis_comm = data.frame(matrix(spec_list_ritis_comm, nrow=attr(spec_list_ritis_comm, "dim")[1]))
      # browser()
      
  #       # missing_spec = unique(c(unlist(spec_list_taxize_sci[spec_list_taxize_sci$MATCH_WORMS!="single match","SEARCHTERM"]), 
  #                  unlist(spec_list_worrms_sci[spec_list_worrms_sci$MATCH_WORMS!="single match","SEARCHTERM"]), 
  #                  unlist(spec_list_ritis_sci[spec_list_ritis_sci$MATCH_ITIS!="single match","SEARCHTERM"])))
  
  
  
  
  #'final assignment of mysterySpec (no matches found) and multimatches by 
  #'removing species that are in the dfinitive list
    mysterySpecCln = mysterySpec[!(mysterySpec$SPEC %in% definitive_Aphia$SPEC) & !(mysterySpec$SPEC %in% definitive_TSN$SPEC),]
    assign("mysterySpec", mysterySpecCln, envir = .GlobalEnv)
    
    
    multiCheck()
}

#     

# 
# 
# }
  #   
  # cat("Check that we're using 'accepted' versions of the ITIS IDs...\n")
  # if (NROW(spec_list_ID[!is.na(spec_list_ID$TSN),])>0){
  #   TSNCheck =  data.frame(submittedtsn=character(),
  #                      acceptedname=character(), 
  #                      acceptedtsn=character(), 
  #                      author=character(), 
  #                    stringsAsFactors=FALSE)
  #  for (i in 1:length(spec_list_ID[!is.na(spec_list_ID$TSN),"TSN"])){
  #    TSNCheck = rbind(TSNCheck,taxize::itis_acceptname(searchtsn = spec_list_ID[!is.na(spec_list_ID$TSN),"TSN"][i]))
  #  } 
  #  spec_list_ID = merge(spec_list_ID, TSNCheck[,c("submittedtsn","acceptedtsn","acceptedname")], by.x = "TSN", by.y = "submittedtsn", all.x=TRUE)
  # }else{
  #   spec_list_ID$acceptedtsn = NA
  #   spec_list_ID$acceptedname = NA
  # }
  # spec_list_ID$TSNFinal = NA
  # spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$TSNFinal <- spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$acceptedtsn
  # spec_list_ID$TSN <- unlist(spec_list_ID$TSNFinal)
  # spec_list_ID$TSNFinal <- NULL
  # spec_list_ID$acceptedtsn <- NULL
  # spec_list_ID$SPEC_SUGGEST <- toupper(spec_list_ID$acceptedname)
  # spec_list_ID$acceptedname <- NULL
  # 
  # # check if we can use any of the tsns we have to get the worms records
  # cat("Trying to find missing WoRMS IDs using found ITIS IDs...\n")
  # for (i in 1:NROW(spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),])){
  #   this <- tryCatch(
  #     {
  #       fromJSON(paste0("http://www.marinespecies.org/rest/AphiaRecordByExternalID/",spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$TSN,"?type=tsn"))$AphiaID
  #     },
  #     error=function(cond){
  #       #message(cond)
  #     }
  #   )
  #   if (!is.null(this)) {
  #     spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$MATCH_WORMS <- "from TSN"
  #     spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$APHIAID<-this
  #   }
  # }
  # 
  # cat("Trying to find missing TSN using found AphiaID...\n")
  # for (k in 1:length(spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"APHIAID"])){
  #   this1 <- tryCatch(
  #     {
  #       worrms::wm_external(as.integer(spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"APHIAID"][k]),type = "tsn")
  #       },
  #     error=function(cond){
  #     }
  #   )
  #   if (!is.null(this1)) {
  #     spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"MATCH_ITIS"][k] <- "from aphiaID"
  #     spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"TSN"]<-this1
  #   }
  # }

#   spec_list_ID[,sci_col]<-toupper(spec_list_ID[,sci_col])
#   spec_list_ID <- merge(df_res, spec_list_ID, by.x = c(sci_col,comm_col), by.y=c(sci_col,comm_col))
#   return(spec_list_ID)
# }
