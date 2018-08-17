#' @title do_taxize
#' @description This function sends either the scientific or common name of each
#' species off to taxize - ids that are returned are then sent to worrms for 
#' verification and alternative sp[elling suggestions]
#' @param df - dataframe with names you want to check against ritis
#' @param chkField - name of the field in the df you want to check
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom taxize get_wormsid
#' @importFrom worrms wm_record
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
do_taxize<-function(df = NULL, 
                     chkField = NULL,
                     logName = logName,
                     searchtype = NULL){
  if (searchtype=="scientific"){
    cat(paste0("\ttaxize > scientific names\n"), file = logName, append = TRUE)
  }else  if (searchtype=="common"){
    cat(paste0("\ttaxize > common names\n"), file = logName, append = TRUE)
  }
  u_df = data.frame(u_rec =unique(df[!is.na(df[chkField]),chkField]),
                    match=NA,
                    multiple_matches = NA,
                    pattern_match = NA,
                    CODE = NA)
  results=df[0,]
  updFields = c("CODE","CODE_SRC","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","SUGG_SPELLING")
  pb <- winProgressBar(title = paste0("APHIAID>TAXIZE>",chkField), label=u_df[1,"u_rec"], min = 0, max = nrow(u_df), width = 300)
  for (i in 1:nrow(u_df)) {
    cat(paste0("\t\ttaxize>",searchtype,">",u_df[i,"u_rec"]), file = logName, append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(u_df[i,"u_rec"]," (", nrow(u_df)-i," left)"))
    this <- tryCatch({
      taxize::get_wormsid(
        query = u_df[i,"u_rec"],
        searchtype = searchtype,
        ask = FALSE,
        verbose = FALSE,
        accepted = FALSE,
        rows = 1,
        messages = FALSE
      )
    },
    error = function(cond) {
    })
    if (is.null(this)){
      cat(paste0(" - NA\n"), file = logName, append = TRUE)
       thisrec = df[df[,chkField]==u_df[i,"u_rec"],]
    }else{
      cat(paste0(" - Results\n"), file = logName, append = TRUE)
      tmp=data.frame(this)
      tmp = tmp[,c("ids","match","multiple_matches","pattern_match")]
      tmp$validID<-NA
      tmp$validSciName<-NA
      tmp$definitive <- FALSE
      for (j in 1:nrow(tmp)){
        thisDefCheck <- tryCatch({
          worrms::wm_record(id=as.numeric(tmp[j,'ids']))
        },
        error = function(cond) {
        })
        if (!is.null(thisDefCheck)& length(thisDefCheck)>0 & !is.null(thisDefCheck$valid_AphiaID) & !is.null(thisDefCheck$valid_name)){
          tmp[j,"validID"]<-thisDefCheck$valid_AphiaID
          tmp[j,"validSciName"]<-thisDefCheck$valid_name
        }
      }  
      if (nrow(tmp)==1){
        tmp$definitive=TRUE
      }
      thisrec = data.frame(match=tmp$match,
                           multiple_matches = tmp$multiple_matches,
                           pattern_match = tmp$pattern_match,
                           u_rec =u_df[i,"u_rec"],
                           CODE = tmp$validID,
                           CODE_SRC = searchtype,
                           CODE_SVC = 'TAXIZE & WORRMS',
                           CODE_TYPE = 'APHIAID',
                           CODE_DEFINITIVE = tmp$definitive,
                           SUGG_SPELLING = trimws(toupper(tmp$validSciName))
      )
      #assess nature of returned code
      if (thisrec$multiple_matches==TRUE & thisrec$pattern_match ==FALSE){
        thisrec$CODE_DEFINITIVE<-TRUE
      }else if (thisrec$multiple_matches==FALSE & thisrec$pattern_match ==FALSE){
        thisrec$CODE_DEFINITIVE<-TRUE
      }else{
        thisrec$CODE_DEFINITIVE<-FALSE
      }
      thisrec$multiple_matches<-thisrec$match<-thisrec$pattern_match<-NULL
      thisrec = merge(df[,-which(colnames(df) %in% updFields)],thisrec, all.y=T, by.x=chkField, by.y = "u_rec")
    }
    results = rbind(results,thisrec)
  }
  close(pb)
  return(results)
}