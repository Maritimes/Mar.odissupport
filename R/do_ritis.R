#' @title do_ritis
#' @description This function sends either the scientific or common name of each
#' species off to ritis
#' @param df - dataframe with names you want to check against ritis
#' @param chkField - name of the field in the df you want to check
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom ritis search_scientific
#' @importFrom ritis usage
#' @importFrom ritis search_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
do_ritis<-function(df = NULL, 
                   chkField = NULL,
                   logName = logName,
                   searchtype = NULL){
  
  cleanTSNs <- function(potents){
    potents$taxonUsageRating<-NA
    potents$accepted_names<-NA
    #if we have a bunch of tsns for a query but don't know which to keep
    for (j in 1:nrow(potents)){
      thisDefCheck <- tryCatch({
        data.frame(ritis::usage(potents[j,"CODE"]))
      },
      error = function(cond) {
      })
      if (!is.null(thisDefCheck)& nrow(thisDefCheck)>0){
        potents[j,"taxonUsageRating"]<-thisDefCheck$taxonUsageRating
        if (nrow(thisDefCheck)>1)browser()
      }
    }  
    
    if (nrow(potents[potents$taxonUsageRating =="valid",])>0)potents = potents[potents$taxonUsageRating =="valid",]
    for (k in 1:nrow(potents)){
      thisDefCheck2 <- tryCatch({
        data.frame(ritis::accepted_names(potents[k,"CODE"]))
      },
      error = function(cond) {
      })
      if (!is.null(thisDefCheck2)& nrow(thisDefCheck2)>0){
        potents[k,"SUGG_SPELLING"]<-trimws(toupper(thisDefCheck2$acceptedName))
        potents[k,"CODE"]<-thisDefCheck2$acceptedTsn
        if (nrow(thisDefCheck2)>1)browser()
      }
    }
    
    if (nrow(potents)==1){
      potents$CODE_DEFINITIVE<-TRUE
      potents$taxonUsageRating<-NULL
      potents$accepted_names<-NULL
      return(potents)
    }  
    invalid_potents <- potents[!(potents$taxonUsageRating %in% "valid"),  ]
    if (nrow(invalid_potents)>0)cat(paste0("opportunity to further limit tsns\n"), file = logName, append = TRUE)
    
    #use taxonusageRating to help
    potents$taxonUsageRating<-NULL
    potents$accepted_names<-NULL
    return(potents)
  }
  if (searchtype=="scientific"){
    cat(paste0("\tritis > scientific names\n"), file = logName, append = TRUE)
  }else  if (searchtype=="common"){
    cat(paste0("\tritis > common names\n"), file = logName, append = TRUE)
  }
  u_df = data.frame(u_rec =unique(df[!is.na(df[chkField]),chkField]),
                    ritisname=NA,
                    CODE = NA)
  results=df[0,]
  updFields = c("CODE","CODE_SRC","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","SUGG_SPELLING")
  pb <- winProgressBar(title = paste0("TSN>ritis>",chkField), label=u_df[1,"u_rec"], min = 0, max = nrow(u_df), width = 300)
  for (i in 1:nrow(u_df)) {
    cat(paste0("\t\tritis>",searchtype,">",u_df[i,"u_rec"]), file = logName, append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(u_df[i,"u_rec"]," (", nrow(u_df)-i," left)"))
    if (searchtype =='scientific'){
      this <- tryCatch({
        as.data.frame(ritis::search_scientific(x = u_df[i,"u_rec"], 
                                               wt = "json", 
                                               raw = FALSE))
      },
      error = function(cond) {
      })
    }else if (searchtype == 'common'){
      this <- tryCatch({
        ritis::search_common(x = u_df[i,"u_rec"], 
                             wt = "json", 
                             raw = FALSE)
      },
      error = function(cond) {
      })
    }
    if (is.null(this)){
      cat(paste0("- NA\n"), file = logName, append = TRUE)
      thisrec = df[df[,chkField]==u_df[i,"u_rec"],]
    }else if (nrow(this)==0){
      cat(paste0("- NA\n"), file = logName, append = TRUE)
      thisrec = df[df[,chkField]==u_df[i,"u_rec"],]
    }else{
      cat(paste0("- Found\n"), file = logName, append = TRUE)
      tmp=unique(data.frame(this))
      #format results so they look the same regardless of searchtype
      if (searchtype=="scientific"){
        tmp=tmp[,c("tsn","combinedName")]
        names(tmp)[names(tmp) == 'combinedName'] <- 'ritisname'
      }
      if (searchtype=="common"){
        tmp=tmp[,c("tsn","commonName")]
        names(tmp)[names(tmp) == 'commonName'] <- 'ritisname'
      }
      thisrec = data.frame(u_rec =u_df[i,"u_rec"],
                           CODE = tmp$tsn,
                           CODE_SRC = searchtype,
                           CODE_SVC = 'RITIS',
                           CODE_TYPE = 'TSN',
                           CODE_DEFINITIVE = FALSE,
                           SUGG_SPELLING = trimws(toupper(tmp$ritisname)))
      
      #assess nature of returned code
      #this prevents genera from returning all known species or a species from returnng a subspecies 
      if (nrow(thisrec[sapply(strsplit(thisrec$SUGG_SPELLING, " "), length) ==sapply(strsplit(u_df[i,"u_rec"], " "), length), ])>0){
        thisrec = thisrec[sapply(strsplit(thisrec$SUGG_SPELLING, " "), length) ==sapply(strsplit(u_df[i,"u_rec"], " "), length), ]  
      }  
      # if (searchtype=="scientific" & nrow(thisrec[thisrec$SUGG_SPELLING==thisrec$u_rec,])==1){
      #   thisrec = thisrec[thisrec$SUGG_SPELLING==thisrec$u_rec,]
      #   thisrec$CODE_DEFINITIVE<-TRUE
      # }
      thisrec = cleanTSNs(thisrec)
      thisrec = merge(df[,-which(colnames(df) %in% updFields)],thisrec, all.y=T, by.x=chkField, by.y = "u_rec")
    }
    results = rbind(results,thisrec)
  }
  close(pb)
  return(results)
}