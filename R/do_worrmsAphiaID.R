#' @title do_worrmsAphiaID
#' @description This function sends TSNS to worrms to see if it can find a 
#' corresponding APHIAID
#' @param recs - vector of taxa names for which we already have a TSN
#' @param knownTSNs - dataframe of taxas for which we have the TSNs and a
#' field called "SCI_COL_CLN" containing and identifying name for the taxa. 
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
#' @importFrom worrms wm_record_by_external
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
do_worrmsAphiaID <-function(recs = NULL, 
                            knownTSNs=NULL,
                        logName = logName){
  total <- length(recs)
  pb <- winProgressBar(title = paste0("APHIAID>WORRMS: via discovered TSNs"), label=recs[1], min = 0, max = total, width = 300)
  
  df= data.frame(joincol = "MMMMMMM",
                 CODE = NA,
                 CODE_SVC = 'WORRMS',
                 CODE_TYPE = 'APHIAID',
                 CODE_SRC = 'TSN',
                 CODE_DEFINITIVE = FALSE,
                 SUGG_SPELLING = NA)
  for (i in 1:total) {
    knownTSNs[knownTSNs$CODE==recs[i],"SCI_COL_CLN"]
    #"(",knownTSNs[knownTSNs$CODE==recs[i],"SCI_COL_CLN"],"/",knownTSNs[knownTSNs$CODE==recs[i],"COMM_COL_CLN"],")
    cat(paste0("\t\t\tworrms>TSN>",recs[i],"\n"), file = logName, append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(knownTSNs[knownTSNs$CODE==recs[i],"SCI_COL_CLN"]," (", total-i," left)"))
    this <- tryCatch(
      {
        worrms::wm_record_by_external(id = as.integer(recs[i]))
        
      },
      error=function(cond){
      }
    )
    if(is.null(this)){
      thisrec= data.frame(joincol = trimws(toupper(recs[i])),
                          CODE = NA,
                          CODE_SVC = 'WORRMS',
                          CODE_TYPE = 'APHIAID',
                          CODE_SRC = 'TSN',
                          CODE_DEFINITIVE = FALSE,
                          SUGG_SPELLING = NA)
    }else{
      thisrec = data.frame(joincol = trimws(toupper(recs[i])),
                           CODE = this$valid_AphiaID,   
                           CODE_SVC = 'WORRMS',
                           CODE_TYPE = 'APHIAID',
                           CODE_SRC = "TSN",
                           CODE_DEFINITIVE = FALSE,
                           SUGG_SPELLING = trimws(toupper(this$valid_name))
      )
      if (length(this$AphiaID)==1)thisrec$CODE_DEFINITIVE<-TRUE
    }
    df = rbind(df,thisrec)
  }
  close(pb)
  df=df[df$joincol!="MMMMMMM",]
  names(df)[names(df) == "joincol"] <- "TSN"
  return(df)
}
