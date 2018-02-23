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
getTaxaIDs <- function(spec_list=NULL, sci_col=NULL, comm_col=NULL){
  #' remove potential confusion from CASE sensitivity
  spec_list = as.data.frame(sapply(spec_list, toupper))

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
  multiMatches <<- data.frame(suggested_name = character(), AphiaID= integer(), AphiaID_src= character(), internal_name = character(), AphiaID_data = character())
  #
  mysterySpec <<-spec_list
  
  if (is.null(comm_col)){
      df_imp<-data.frame("tmpzzz" = spec_list[,sci_col])
      names(df_imp)[names(df_imp)=="tmpzzz"]<-sci_col
  } else {
    df_imp<-data.frame(spec_list[,c(sci_col, comm_col)])
  }

  df_res<-spec_list

  # #function for modifying case of entries
  # proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), perl = TRUE)
  # 
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

# Compare results of services that look for the same ID -------------------
  matchCheck <-function(spec_list, df1, df2, searchtype = NULL, IDtype = NULL){
    if (IDtype =="AphiaID"){
      theID = "AphiaID"
      theIDComment = "AphiaID_src"
    }
    if (searchtype =="scientific"){
      searchtypeVal = "scientific name"
      colname = "SPEC"
    }else{
      searchtypeVal = "common name"
      colname = "COMM"
    }
    #'the block below exists with the assumption that the only 2 groups of 
    #'results that will be sent to this function are taxize and worrms.
    #'This assumption and the block may need to change in the future.
    if(any(grep("taxize",colnames(df1)))){
      df1_ID = "ID_taxize"
      df1_MATCH = "MATCH_taxize"
      df1_type = "Taxize"
      df2_ID = "ID_worrms"
      df2_MATCH = "MATCH_worrms"
      df2_type = "Worrms"
    }else if (any(grep("worrms",colnames(df1)))){
      df1_ID = "ID_worrms"
      df1_MATCH = "MATCH_worrms"
      df1_type = "Worrms"
      df2_ID = "ID_taxize"
      df2_MATCH = "MATCH_taxize"
      df2_type = "Taxize"
  }
    keeperStrings = c("single match","only accepted")
    results = merge(df1, df2, by="SEARCHTERM")
    results[theID] = NA
    results[theIDComment] = NA
    results[grep(paste(keeperStrings,collapse="|"),results[,df1_MATCH]),theID]<-results[grep(paste(keeperStrings,collapse="|"),results[,df1_MATCH]),df1_ID]
    results[grep(paste(keeperStrings,collapse="|"),results[,df2_MATCH]),theID]<-results[grep(paste(keeperStrings,collapse="|"),results[,df2_MATCH]),df2_ID]
    #no ID from either df = unmatched
    if(nrow(results[is.na(results[,df1_ID]) & (is.na(results[,df2_ID])),])>0)
      results[is.na(results[,df1_ID]) & (is.na(results[,df2_ID])),theIDComment]<-"Unmatched"
    #only 1 id present, and it was in the keeper strings (i.e. a confident match)
    if(nrow(results[(!is.na(results[,df1_ID]) & is.na(results[,df2_ID])) & !is.na(results[,theID]),])>0)
      results[(!is.na(results[,df1_ID]) & is.na(results[,df2_ID])) & !is.na(results[,theID]),theIDComment] <-paste(df1_type, "Only")
    if(nrow(results[(is.na(results[,df1_ID]) & !is.na(results[,df2_ID])) & !is.na(results[,theID]),])>0)
      results[(is.na(results[,df1_ID]) & !is.na(results[,df2_ID])) & !is.na(results[,theID]),theIDComment] <-paste(df2_type, "Only")
    #both ids present, and at least one was in the keeper strings (i.e. a confident match)
    if(nrow(results[(!is.na(results[,df1_ID]) & (!is.na(results[,df2_ID])) & (!is.na(results[,theID]))) & results[,df1_ID] ==results[,df2_ID],])>0)
      results[(!is.na(results[,df1_ID]) & (!is.na(results[,df2_ID])) & (!is.na(results[,theID]))) & results[,df1_ID] ==results[,df2_ID],theIDComment] <-paste(df1_type, "&",df2_type)
    #all that remain should be those which had 1 or more matches, but none within the keeper strings
    if(nrow(results[is.na(results[,theIDComment]),])>0)
      results[is.na(results[,theIDComment]),]<-"Ambiguous Match"
    
    results["AphiaID_data"] =NA
    if (nrow(results[!is.na(results$AphiaID),])>0)
      results[!is.na(results$AphiaID),"AphiaID_data"]<-searchtypeVal
    results <-results[,c("SEARCHTERM",theID,theIDComment,"AphiaID_data")]
    names(results)[names(results) == 'SEARCHTERM'] <- colname
    #get the records we have confident IDs for
    thesedefinitive_Aphia = results[!(results$AphiaID_src %in% c("Ambiguous Match", "Unmatched")),]
    #drop any that are already in our definitve list

    thesedefinitive_Aphia = thesedefinitive_Aphia[!(thesedefinitive_Aphia$SPEC %in% definitive_Aphia$SPEC),]
    definitive_Aphia<<-rbind(definitive_Aphia,thesedefinitive_Aphia)
    #remove recs from multimatches that have been accepted
       # if (nrow(definitive_Aphia)>0){
          moreChecks <-merge(spec_list[!(spec_list[,colname] %in% definitive_Aphia$SPEC),], results[,c(colname,"AphiaID","AphiaID_src")])#,by.x=colname, by.y = "SPEC")
       #  }
    # else{
    #      moreChecks <- merge(spec_list[!(spec_list[,colname] %in% results[,colname]),], results[,c("SPEC","AphiaID","AphiaID_src")],by.x=colname, by.y = "SPEC")
    #     }
    return(moreChecks)
  }

# Tests against various services ------------------------------------------
  checkAndCompareAphia<-function(df= NULL, searchtype = NULL, IDtype = "AphiaID"){
    # check for matches based on sci name, then common name (if provided)
    cat(paste("Checking using ",searchtype," names\n"))
    cat("\tagainst Taxize (AphiaID)....\n")
    spec_list_taxize = as.data.frame(t(sapply_pb(df,chk_taxize)))
    spec_list_taxize = data.frame("SEARCHTERM" = toupper(unlist(spec_list_taxize$SEARCHTERM)),
                                  "ID_taxize" = unlist(spec_list_taxize$ID_taxize),
                                  "MATCH_taxize" = unlist(spec_list_taxize$MATCH_taxize))
    cat("\tagainst worrms (AphiaID)\n")
    spec_list_worrms = as.data.frame(t(sapply_pb(df,chk_worrms)))
    spec_list_worrms = data.frame("SEARCHTERM" = toupper(unlist(spec_list_worrms$SEARCHTERM)),
                                  "ID_worrms" = unlist(spec_list_worrms$ID_worrms),
                                  "MATCH_worrms" = unlist(spec_list_worrms$MATCH_worrms))
    
    #combine all of the results from the scientific names
    #lets deal with aphiaids first
    doMore = matchCheck(spec_list = spec_list,
                        df1 = spec_list_taxize, df2 = spec_list_worrms,
                        searchtype = searchtype, IDtype = IDtype)
    if (nrow(doMore)==0)cat("\nAll IDs found\n")
    return(doMore)
  }
# Functions for checking web services -------------------------------------
  chk_taxize <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use taxize functions to check stuff
    this <- tryCatch(
      {
        taxize::get_wormsid(query=searchterm, searchtype=searchtype, ask=ask, verbose=verbose, accepted = TRUE,rows = 1)
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
    return(this)
  }

  chk_ritis <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use ritis functions to check stuff
    this = resultFormat
    for (i in 1:NROW(searchterm)){
      if (searchtype == 'scientific'){
        resName = "combinedName"
        thisi =  as.data.frame(ritis::search_scientific(searchterm[i], wt = "json", raw = FALSE))
      }else{
        resName = "commonName"
        thisi = as.data.frame(ritis::search_common(searchterm[i], wt = "json", raw = FALSE))
      }
      
      if (nrow(thisi)==0) {
        thisi=data.frame("SEARCHTERM" = searchterm[i], "ID_ritis" = NA, "MATCH_ritis" = "not found")
      }else if (nrow(thisi)==1) {
        thisi = cbind("SEARCHTERM"= searchterm[i],"ID_ritis" = thisi$tsn, "MATCH_ritis" = "single match")
      }else{
        mult = thisi[,c(resName,"TSN")]
        colnames(mult)<-c("suggested_name","TSN")
        mult$TSN_internal_name = searchterm[i]
        mult$TSN_src = "ritis"
        mult$TSN_data = searchtype
        multiMatches <<- rbind(multiMatches, mult)
          if (length(thisi[thisi[,resName] == searchterm[i],"tsn"])==1){
            thisi = cbind(searchterm[i], thisi[thisi[,resName] == searchterm[i],"tsn"], "multi match - retained identical spelling")
          }else{
            thisi = cbind(searchterm[i], thisi$tsn[1], "multi match - no identical spelling; retained first")
          }
      }

      this = rbind(this, thisi)
    }
    this = this[rowSums(is.na(this)) != ncol(this),]
    return(this)
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
    chk_worrms <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    #use worrms functions to check stuff
    this = resultFormat
    for (i in 1:NROW(searchterm)){
      resName = "scientificname"
      if (searchtype == 'scientific'){
        thisi <- tryCatch(
          {
            worrms::wm_records_name(searchterm[i], marine_only = T)
          },
          error=function(cond){
          }
        )
      }else{
        thisi <- tryCatch(
          {
            worrms::wm_records_common(searchterm[i], marine_only = T)
          },
          error=function(cond){
            #message(cond)
          }
        )
      }
      crapRec = data.frame("SEARCHTERM" = searchterm[i], "ID_worrms" = NA, "MATCH_worrms" = "not found")
      if (is.null(thisi)) {
        thisi <- crapRec
      } else if (nrow(thisi)==0) {
        thisi <- crapRec
      }else if (nrow(thisi)==1) {
        thisi = cbind("SEARCHTERM"= searchterm[i],"ID_worrms" = thisi$AphiaID, "MATCH_worrms" = "single match")
      }else{
        mult = thisi[,c(resName,"AphiaID")]
        colnames(mult)<-c("suggested_name","AphiaID")
        mult$AphiaID_src = "worrms"
        mult$AphiaID_internal_name = searchterm[i]
        mult$AphiaID_data = searchtype
        multiMatches <<- rbind(multiMatches, mult)
        
        if (nrow(thisi[thisi$status=="accepted",])==1){
          thisi = cbind("SEARCHTERM" = searchterm[i],"ID_worrms" = thisi[thisi$status=="accepted","AphiaID"], "MATCH_worrms" = "multi match - kept only accepted")
          colnames(thisi)[colnames(thisi)=="AphiaID"] <- "ID_worrms"
        } else if (nrow(thisi[thisi$status=="accepted",])>1) {
          thisi = cbind("SEARCHTERM" = searchterm[i],"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - kept 1st accepted")
        }else {
          #multiMatches =rbind(multiMatches, data.frame(sci = thisi$scientificname, src = "worrms", ID = as.integer(thisi$AphiaID)))
          thisi = cbind("SEARCHTERM" = searchterm[i],"ID_worrms" = thisi$AphiaID[1], "MATCH_worrms" = "multi match - none accepted; retained first")
        }

      }
      this = rbind(this, thisi)
    }
    return(this)
  }

# Start of actual check (all functions above) -----------------------------
  #Check for aphias by scientific name
  sci_missing = checkAndCompareAphia(df = df_imp[,sci_col],searchtype="scientific")
  # browser()
  #Check for aphias by common name  
  if (nrow(sci_missing)>0){
    cat("\nChecking ambiguously and unmatched records using Common names\n")
    comm_results = checkAndCompareAphia(df = sci_missing[,comm_col],searchtype="common")
  }
  
 
  cat("\nFinding TSNs using the 'definitive' AphiaIDs\n")
  TSNsFromAphias = sapply_pb(definitive_Aphia$AphiaID,chk_worrmsTSN)
  TSNsFromAphias = do.call(rbind,TSNsFromAphias)
  TSNsFromAphias = data.frame(join=rownames(TSNsFromAphias), data=TSNsFromAphias, row.names = NULL)
  names(TSNsFromAphias)<-c("AphiaID", "TSN")
  definitive_TSN = data.frame(SPEC = definitive_Aphia$SPEC, TSN_src = "worrms", TSN_data = "AphiaId", AphiaID = definitive_Aphia$AphiaID, row.names = NULL)
  definitive_TSN = merge(definitive_TSN, TSNsFromAphias, all.x = T)
  

  multitest = sapply_pb(multiMatches$AphiaID,chk_worrmsTSN)
  multitest = do.call(rbind,multitest)
  multitest = data.frame(join=rownames(multitest), data=multitest, TSN_src = "worrms", TSN_data ="AphiaId", row.names = NULL)
  names(multitest)<-c("AphiaID", "TSN","TSN_src","TSN_data")
  multiMatches = merge(multiMatches,multitest, all.x = T, by.x="AphiaID", by.y = "AphiaID")
  browser() 
  #find tsn using ritis
  # definitive_Aphia$TSN<-sapply_pb()
  # browser()l

  # definitive $TSN <- worrms::wm_external_(id = as.integer(definitive$AphiaID))
  # definitive[!is.na(definitive$TSN),"TSN_src"]<-"Worrms via AphiaID"
  
  #'final assignment of mysterySpec (no matches found) and multimatches by 
  #'removing species that are in the dfinitive list

    mysterySpecCln = mysterySpec[!(mysterySpec$SPEC %in% definitive_Aphia$SPEC),]
    assign("mysterySpec", mysterySpecCln, envir = .GlobalEnv)
    clnMulti = multiMatches[!(toupper(multiMatches$AphiaID_internal_name) %in% definitive_Aphia$SPEC),]
    assign("multiMatches", clnMulti, envir = .GlobalEnv)
    browser()
}

#     
#     cat("\tagainst ritis (TSN) \n")
#     spec_list_ritis_sci = as.data.frame(t(sapply_pb(df_imp[,sci_col],chk_ritis)))
#     spec_list_ritis_sci = data.frame("SEARCHTERM" = toupper(unlist(spec_list_ritis_sci$V1)), 
#                                      "ID_ritis" = unlist(spec_list_ritis_sci$V2), 
#                                      "MATCH_ritis" = unlist(spec_list_ritis_sci$V3))
# }
#       # missing_spec = unique(c(unlist(spec_list_taxize_sci[spec_list_taxize_sci$MATCH_WORMS!="single match","SEARCHTERM"]), 
  #                  unlist(spec_list_worrms_sci[spec_list_worrms_sci$MATCH_WORMS!="single match","SEARCHTERM"]), 
  #                  unlist(spec_list_ritis_sci[spec_list_ritis_sci$MATCH_ITIS!="single match","SEARCHTERM"])))
    

#     cat("\tagainst ritis (TSN) \n")
#     spec_list_ritis_comm = as.data.frame(t(sapply_pb(df_imp[,comm_col],chk_ritis, searchtype="common")))
#     spec_list_ritis_comm = data.frame(matrix(spec_list_ritis_comm, nrow=attr(spec_list_ritis_comm, "dim")[1]))
# browser()
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
