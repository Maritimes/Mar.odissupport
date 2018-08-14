#' @title do_worrmsTSN
#' @description This function sends APHIAIDs to worrms to see if it can find a 
#' corresponding TSN
#' @param recs - vector of taxa names for which we already have an APHIAID
#' @param knownAphias - dataframe of taxas for which we have the AphiaIDs and a
#' field called "SCI_COL_CLN" containing and identifying name for the taxa. 
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
#' @importFrom worrms wm_external
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
do_worrmsTSN <-function(recs = NULL, 
                         knownAphias=NULL,
                         logName = logName){
  total <- length(recs)
  pb <- winProgressBar(title = paste0("TSN>WORRMS: via discovered APHIAIDs"), label=recs[1], min = 0, max = total, width = 300)
  
  df= data.frame(joincol = "MMMMMMM",
                 CODE = NA,
                 CODE_SVC = 'WORRMS',
                 CODE_TYPE = 'TSN',
                 CODE_SRC = 'APHIAID',
                 CODE_DEFINITIVE = FALSE,
                 SUGG_SPELLING = NA)
  # names(df)[names(df)=="CODE"]= "AphiaID"
  # names(df)[names(df)=="CODE_DEFINITIVE"]= "AphiaID_DEF"
  
  #added loop since one record with no results botched the whole call;
  #likely a huge performance hit

  for (i in 1:total) {
    knownAphias[knownAphias$CODE==recs[i],"SCI_COL_CLN"]
    cat(paste0("\t\tworrms>APHIAID>",recs[i],"(",knownAphias[knownAphias$CODE==recs[i],"SCI_COL_CLN"],"/",knownAphias[knownAphias$CODE==recs[i],"COMM_COL_CLN"],")\n"), file = logName, append = TRUE)
      setWinProgressBar(pb, i, title = NULL, label = paste0(knownAphias[knownAphias$CODE==recs[i],"SCI_COL_CLN"]," (", total-i," left)"))
    this <- tryCatch(
      {
        worrms::wm_external(id = as.integer(recs[i]))
      },
      error=function(cond){
      }
    )
    if(is.null(this)){
      thisrec= data.frame(joincol = trimws(toupper(recs[i])),
                       CODE = NA,
                       CODE_SVC = 'WORRMS',
                       CODE_TYPE = 'TSN',
                       CODE_SRC = 'APHIAID',
                       CODE_DEFINITIVE = FALSE,
                       SUGG_SPELLING = NA)
    }else{
      
      thisrec = data.frame(joincol = trimws(toupper(recs[i])),
                        CODE = this,   
                        CODE_SVC = 'WORRMS',
                        CODE_TYPE = 'TSN',
                        CODE_SRC = "APHIAID",
                        CODE_DEFINITIVE = FALSE,
                        SUGG_SPELLING = NA
      )
      for (j in 1:nrow(thisrec)){
        thisDefCheck <- tryCatch({
          data.frame(ritis::usage(thisrec[j,"CODE"]))
        },
        error = function(cond) {
        })
        if(is.null(thisDefCheck)){
          thisrec[j,"CODE_DEFINITIVE"]<-FALSE
        }else{
          thisrec[j,"CODE_DEFINITIVE"]<-ifelse(thisDefCheck["taxonUsageRating"]=="valid",TRUE,FALSE)
        }
      }
      
    }
    df = rbind(df,thisrec)
  }
  close(pb)
  df=df[df$joincol!="MMMMMMM",]
  names(df)[names(df) == "joincol"] <- "APHIAID"
  return(df)
}
