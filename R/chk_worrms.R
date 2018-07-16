#' @title chk_worrms
#' @description This function sends either the scientific or common name of each
#' species off to worrms
#' @param recs -  - vector of taxa names for which the APHIAID is either unknown or not
#' definitive
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom worrms wm_records_name
#' @importFrom worrms wm_records_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_worrms <- function(recs = NULL,
                       searchtype = NULL) {
  col = ifelse(searchtype=="scientific","SCI_COL_CLN", "COMM_COL_CLN")
  df= data.frame(joincol = "MMMMMMM",
                   CODE = NA,
                   CODE_SVC = 'WORRMS',
                   CODE_TYPE = 'APHIAID',
                   CODE_SRC = searchtype,
                   CODE_DEFINITIVE = FALSE,
                   SUGG_SPELLING = NA)
  total <- length(recs)

  pb <- winProgressBar(title = paste0("APHIAID>WORRMS: via ", searchtype," names"), label=recs[1], min = 0, max = total, width = 300)
  df=df[0,]
  for (i in 1:total) {
    cat(paste0("\tworrms|",searchtype,"|",recs[i],"\n"), file = "getTaxaIDs.log", append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(recs[i]," (", total-i," left)"))
    if (searchtype == 'scientific') {
      this <- tryCatch({
        worrms::wm_records_name(recs[i], fuzzy = F, marine_only = T)
      },
      error = function(cond) {
      })
    } else if (searchtype == 'common'){
      this <- tryCatch({
        worrms::wm_records_common(recs[i], fuzzy = F, marine_only = T)
      },
      error = function(cond) {
      })
    }
    if(is.null(this)){
      thisrec= data.frame(joincol = trimws(toupper(recs[i])),
                       CODE = NA,
                       CODE_SVC = 'WORRMS',
                       CODE_TYPE = 'APHIAID',
                       CODE_SRC = searchtype,
                       CODE_DEFINITIVE = FALSE,
                       SUGG_SPELLING = NA)
    } else {
      tmp = data.frame(this[,c("status","valid_AphiaID","valid_name")]) #"unacceptreason"
      #if (length(tmp[tmp$status %in% c("unaccepted") & !is.na(tmp$valid_AphiaID),"valid_AphiaID"]>0))browser()
      thisrec = data.frame(
                        status = tmp$status,
                        joincol = trimws(toupper(recs[i])),
                        CODE = tmp$valid_AphiaID,
                        CODE_SRC = searchtype,
                        CODE_SVC = 'WORRMS',
                        CODE_TYPE = 'APHIAID',
                        CODE_DEFINITIVE = FALSE,
                        SUGG_SPELLING =trimws(toupper(tmp$valid_name))
      )
      thisrec[thisrec$status %in% c("accepted") & !is.na(thisrec$CODE),"CODE_DEFINITIVE"]<-TRUE 
      #unaccepted name, but valid code and suggested spelling
      thisrec[thisrec$status %in% c("unaccepted") & !is.na(thisrec$CODE) & !is.na(thisrec$SUGG_SPELLING),"CODE_DEFINITIVE"]<-TRUE
      
      thisrec$status<-NULL
      rm(tmp)
    }   
      df = rbind(df,thisrec)
  }
  close(pb)
  df=df[df$joincol!="MMMMMMM",]
  names(df)[names(df) == "joincol"] <- col
  return(df)
}
