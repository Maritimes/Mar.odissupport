#' @title do_worrms
#' @description This function sends either the scientific or common name of each
#' species off to worrms
#' @param df - dataframe with names you want to check against ritis
#' @param chkField - name of the field in the df you want to check
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom worrms wm_records_name
#' @importFrom worrms wm_records_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
do_worrms<-function(df = NULL, 
                    chkField = NULL,
                    logName = logName,
                    searchtype = NULL){
  if (searchtype=="scientific"){
    cat(paste0("\tworrms > scientific names\n"), file = logName, append = TRUE)
  }else  if (searchtype=="common"){
    cat(paste0("\tworrms > common names\n"), file = logName, append = TRUE)
  }
  u_df = data.frame(u_rec =unique(df[!is.na(df[chkField]),chkField]),
                    valid_AphiaID = NA,
                    valid_name = NA,
                    CODE = NA)
  results=df[0,]
  updFields = c("CODE","CODE_SRC","CODE_SVC","CODE_TYPE","CODE_DEFINITIVE","SUGG_SPELLING")
  pb <- winProgressBar(title = paste0("APHIAID>WORRMS>",chkField), label=u_df[1,"u_rec"], min = 0, max = nrow(u_df), width = 300)
  for (i in 1:nrow(u_df)) {
    cat(paste0("\t\tworrms>",searchtype,">",u_df[i,"u_rec"]), file = logName, append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(u_df[i,"u_rec"]," (", nrow(u_df)-i," left)"))
    if (searchtype == 'scientific') {
      this <- tryCatch({
        worrms::wm_records_name(
          name = u_df[i,"u_rec"],
          fuzzy = F, 
          marine_only = F
        )
      },
      error = function(cond) {
      })
    } else if (searchtype == 'common'){
      this <- tryCatch({
        worrms::wm_records_common(u_df[i,"u_rec"], 
                                  fuzzy = F, 
                                  marine_only = F
        )
      },
      error = function(cond) {
      })
    }
    if (is.null(this)){
      cat(paste0(" - NA\n"), file = logName, append = TRUE)
      thisrec = df[df[,chkField]==u_df[i,"u_rec"],]
    }else{
      cat(paste0(" - Results\n"), file = logName, append = TRUE)
      tmp=data.frame(this)
      tmp = tmp[,c("status","valid_AphiaID","valid_name","scientificname")]
      thisrec = data.frame(status=tmp$status,
                           u_rec =u_df[i,"u_rec"],
                           CODE = tmp$valid_AphiaID,
                           CODE_SRC = searchtype,
                           CODE_SVC = 'WORRMS',
                           CODE_TYPE = 'APHIAID',
                           CODE_DEFINITIVE = FALSE,
                           SUGG_SPELLING = trimws(toupper(tmp$valid_name)))
      # #assess nature of returned code
      #if one record is returned with a valid code that's pretty definitive 
      #status like "unaccepted" are still accompanied by code and spelling of the valid alternative
      if (nrow(thisrec[!is.na(thisrec$CODE),])==1) thisrec = thisrec[!is.na(thisrec$CODE),]
      if (nrow(thisrec) ==1) {
        thisrec$CODE_DEFINITIVE <- TRUE
      }else if (nrow(thisrec) >1){
        if (nrow(thisrec[thisrec$status == "accepted", ])>0){
          thisrec = thisrec[thisrec$status == "accepted", ]
          if (nrow(thisrec)==1)thisrec$CODE_DEFINITIVE <- TRUE
        }
        if (searchtype=="common" & "SCI_COL_CLN" %in% names(df)){
          thisSpec = df[which(df$COMM_COL_CLN ==u_df[i,"u_rec"]),"SCI_COL_CLN"]
          #if scientific name exists and matches one record with what the common
          #name returned, this is definitive
          #if there are several records that match, limit to these
          if (!is.na(thisSpec)){
            if (nrow(thisrec[thisrec$SUGG_SPELLING==thisSpec,])==1){
              cat("Multiple returns resolved\n", file = logName, append = TRUE)
              thisrec$SUGG_SPELLING<-thisSpec
              thisrec$CODE_DEFINITIVE <- TRUE
            }else if (nrow(thisrec[thisrec$SUGG_SPELLING==thisSpec,])>1){
              cat("Multiple returns for this\n", file = logName, append = TRUE)
              thisrec = thisrec[thisrec$SUGG_SPELLING==thisSpec,]
            }
          }
        }
      }
      thisrec$status<-NULL
      thisrec = merge(df[,-which(colnames(df) %in% updFields)],thisrec, all.y=T, by.x=chkField, by.y = "u_rec")
    }
    results = rbind(results,thisrec)
  }
  close(pb)
  return(results)
}