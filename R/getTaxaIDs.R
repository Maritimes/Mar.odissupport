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
  if (is.null(comm_col)){
      df_imp<-data.frame("tmpzzz" = spec_list[,sci_col])
      names(df_imp)[names(df_imp)=="tmpzzz"]<-sci_col
  } else {
    df_imp<-data.frame(spec_list[,c(sci_col, comm_col)])
  }

  df_res<-spec_list

# Functions ---------------------------------------------------------------
  sapply_pb <- function(X, FUN, ...)
  {
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

  proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), perl = TRUE)

  chk_WoRMS <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    this = taxize::get_wormsid(query=searchterm, searchtype=searchtype, ask=ask, verbose=verbose)
    this = as.data.frame(cbind(SEARCHTERM = searchterm, APHIAID=as.vector(this),MATCH_WORMS=attr(this,"match"), METHOD_WORMS = searchtype))
    if (NROW(this[grep("NA due to ask=FALSE", this$MATCH_WORMS),])>0)  this[grep("NA due to ask=FALSE", this$MATCH_WORMS),]$MATCH_WORMS<-"multi match"
    return(this)
  }

  chk_ITIS <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    this = data.frame(SEARCHTERM = NA, TSN = NA, MATCH_ITIS = NA, METHOD_ITIS = NA)
    for (i in 1:NROW(searchterm)){
      if (searchtype == 'scientific'){
        thisi = ritis::search_scientific(searchterm[i], wt = "json", raw = FALSE)
      }else{
        thisi = ritis::search_common(searchterm[i], wt = "json", raw = FALSE)
      }
      if (nrow(thisi)==0) {
        thisi=data.frame(SEARCHTERM = searchterm[i], TSN = NA, MATCH_ITIS = "not found", METHOD_ITIS = NA)
      }else if (nrow(thisi)==1) {
        thisi = cbind(SEARCHTERM= searchterm[i],TSN = thisi$tsn, MATCH_ITIS = "found", METHOD_ITIS = searchtype)
      }else{
        thisi = cbind(SEARCHTERM = searchterm[i],TSN = NA, MATCH_ITIS = "multi match", METHOD_ITIS = searchtype)
        #could also populate field with vector of matches, but this takes it out
        #of sync with the worms data and some spp get MANY matches (e.g. 100s),
        #TSN = paste0(thisi$tsn, collapse=", ")
      }
      this = rbind(this, thisi)
    }
    this = this[rowSums(is.na(this)) != ncol(this),]
    return(this)
  }

# Processing Steps --------------------------------------------------------
    df_imp[,sci_col] = proper(df_imp[,sci_col])
  # WoRMS
  # check for matches based on sci name, then common name (if provided)
  cat("Checking WoRMS for Scientific names\n")
  chk_WoRMSspec_list_sci = as.data.frame(t(sapply_pb(df_imp[,sci_col],chk_WoRMS)))
  chk_WoRMSspec_list_sci_miss = chk_WoRMSspec_list_sci[is.na(chk_WoRMSspec_list_sci$APHIAID),]

  if (NROW(chk_WoRMSspec_list_sci_miss)>0){
    if (!is.null(comm_col)) {
      cat("Checking WoRMS for common names\n")
      chk_WoRMSspec_list_comm = as.data.frame(t(sapply_pb(chk_WoRMSspec_list_sci_miss$SEARCHTERM, chk_WoRMS, searchtype="common")))
      chk_WoRMSspec_list_comm = merge(df_imp, chk_WoRMSspec_list_comm, by.x = comm_col, by.y = "SEARCHTERM")
      spec_list_WoRMS = merge(df_imp, chk_WoRMSspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
      spec_list_WoRMS = rbind(spec_list_WoRMS, chk_WoRMSspec_list_comm)
      spec_list_WoRMS = data.frame("tmp1" = spec_list_WoRMS[,sci_col],
                                   "tmp2" = spec_list_WoRMS[,comm_col],
                                   "APHIAID"=unlist(spec_list_WoRMS$APHIAID),
                                   "MATCH_WORMS"=unlist(spec_list_WoRMS$MATCH_WORMS ),
                                   "METHOD_WORMS"=unlist(spec_list_WoRMS$METHOD_WORMS))
      names(spec_list_WoRMS)[names(spec_list_WoRMS) == "tmp1"] <- sci_col
      names(spec_list_WoRMS)[names(spec_list_WoRMS) == "tmp2"] <- comm_col
    }else{
      spec_list_WoRMS = merge(df_imp, chk_WoRMSspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
      spec_list_WoRMS = data.frame("tmp1" = spec_list_WoRMS[,sci_col],
                                   "APHIAID"=unlist(spec_list_WoRMS$APHIAID),
                                   "MATCH_WORMS"=unlist(spec_list_WoRMS$MATCH_WORMS ),
                                   "METHOD_WORMS"=unlist(spec_list_WoRMS$METHOD_WORMS))
      names(spec_list_WoRMS)[names(spec_list_WoRMS) == "tmp1"] <- sci_col
    }
  }else{
    spec_list_WoRMS = merge(df_imp, chk_WoRMSspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
  }

  # ITIS
  # check for matches based on sci name, then common name (if provided)
  cat("Checking ITIS for Scientific names\n")
  chk_ITISspec_list_sci =  data.frame(SEARCHTERM=character(),
                          TSN=integer(), 
                          MATCH_ITIS=character(), 
                          METHOD_ITIS=character(), 
                         stringsAsFactors=FALSE)
  for (j in 1:length(df_imp[,sci_col])){
    chk_ITISspec_list_sci = rbind(chk_ITISspec_list_sci, chk_ITIS(df_imp[,sci_col][j]))
  }
  
  # chk_ITISspec_list_sci = as.data.frame(t(sapply_pb(df_imp[,sci_col],chk_ITIS)))
  chk_ITISspec_list_sci_miss = chk_ITISspec_list_sci[is.na(chk_ITISspec_list_sci$TSN),]

  if (NROW(chk_ITISspec_list_sci_miss)>0){

    if (!is.null(comm_col)) {
      cat("Checking ITIS for common names\n")
      chk_ITISspec_list_comm = as.data.frame(t(sapply_pb(chk_ITISspec_list_sci_miss$SEARCHTERM, chk_ITIS, searchtype="common")))
      chk_ITISspec_list_comm = merge(df_imp, chk_ITISspec_list_comm, by.x = comm_col, by.y = "SEARCHTERM")
      spec_list_ITIS = merge(df_imp, chk_ITISspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
      spec_list_ITIS = rbind(spec_list_ITIS, chk_ITISspec_list_comm)
      spec_list_ITIS = data.frame("tmp1"=spec_list_ITIS[,sci_col],
                                  "tmp2"=spec_list_ITIS[,comm_col],
                                  "TSN"=unlist(spec_list_ITIS$TSN),
                                  "MATCH_ITIS"=unlist(spec_list_ITIS$MATCH_ITIS ),
                                  "METHOD_ITIS"=unlist(spec_list_ITIS$METHOD_ITIS))
      names(spec_list_ITIS)[names(spec_list_ITIS) == "tmp1"] <- sci_col
      names(spec_list_ITIS)[names(spec_list_ITIS) == "tmp2"] <- comm_col
  }else{
    spec_list_ITIS = merge(df_imp, chk_ITISspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
    spec_list_ITIS = data.frame("tmp1"=spec_list_ITIS[,sci_col],
                                "TSN"=unlist(spec_list_ITIS$TSN),
                                "MATCH_ITIS"=unlist(spec_list_ITIS$MATCH_ITIS ),
                                "METHOD_ITIS"=unlist(spec_list_ITIS$METHOD_ITIS))
    names(spec_list_ITIS)[names(spec_list_ITIS) == "tmp1"] <- sci_col
  }
  }else{
    spec_list_ITIS = merge(df_imp, chk_ITISspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM", all.x=TRUE)
  }
  spec_list_ID <- merge(spec_list_WoRMS, spec_list_ITIS, by=c(sci_col, comm_col))

  cat("Check that we're using 'accepted' versions of the ITIS IDs...\n")
  if (NROW(spec_list_ID[!is.na(spec_list_ID$TSN),])>0){
    TSNCheck =  data.frame(submittedtsn=character(),
                       acceptedname=character(), 
                       acceptedtsn=character(), 
                       author=character(), 
                     stringsAsFactors=FALSE)
   for (i in 1:length(spec_list_ID[!is.na(spec_list_ID$TSN),"TSN"])){
     TSNCheck = rbind(TSNCheck,taxize::itis_acceptname(searchtsn = spec_list_ID[!is.na(spec_list_ID$TSN),"TSN"][i]))
   } 
   spec_list_ID = merge(spec_list_ID, TSNCheck[,c("submittedtsn","acceptedtsn","acceptedname")], by.x = "TSN", by.y = "submittedtsn", all.x=TRUE)
  }else{
    spec_list_ID$acceptedtsn = NA
    spec_list_ID$acceptedname = NA
  }
  spec_list_ID$TSNFinal = NA
  spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$TSNFinal <- spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$acceptedtsn
  spec_list_ID$TSN <- unlist(spec_list_ID$TSNFinal)
  spec_list_ID$TSNFinal <- NULL
  spec_list_ID$acceptedtsn <- NULL
  spec_list_ID$SPEC_SUGGEST <- toupper(spec_list_ID$acceptedname)
  spec_list_ID$acceptedname <- NULL

  # check if we can use any of the tsns we have to get the worms records
  cat("Trying to find missing WoRMS IDs using found ITIS IDs...\n")
  for (i in 1:NROW(spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),])){
    this <- tryCatch(
      {
        fromJSON(paste0("http://www.marinespecies.org/rest/AphiaRecordByExternalID/",spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$TSN,"?type=tsn"))$AphiaID
      },
      error=function(cond){
        #message(cond)
      }
    )
    if (!is.null(this)) {
      spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$MATCH_WORMS <- "from TSN"
      spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$APHIAID<-this
    }
  }

  cat("Trying to find missing TSN using found AphiaID...\n")
  for (k in 1:length(spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"APHIAID"])){
    this1 <- tryCatch(
      {
        wm_external(as.integer(spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"APHIAID"][k]),type = "tsn")
        },
      error=function(cond){
      }
    )
    if (!is.null(this1)) {
      spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"MATCH_ITIS"][k] <- "from aphiaID"
      spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),"TSN"]<-this1
    }
  }

  spec_list_ID[,sci_col]<-toupper(spec_list_ID[,sci_col])
  spec_list_ID <- merge(df_res, spec_list_ID, by.x = c(sci_col,comm_col), by.y=c(sci_col,comm_col))
  return(spec_list_ID)
}
