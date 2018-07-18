#' @title chk_taxize
#' @description This function sends either the scientific or common name of each
#' species off to taxize
#' @param recs - - vector of taxa names for which the APHIAID is either unknown or not
#' definitive
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom taxize get_wormsid
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_taxize <- function(recs = NULL,
           searchtype = NULL) {
  col = ifelse(searchtype=="scientific","SCI_COL_CLN", "COMM_COL_CLN")
  df= data.frame(joincol = "MMMMMMM",
                   CODE = NA,
                   CODE_SVC = 'TAXIZE',
                   CODE_TYPE = 'APHIAID',
                   CODE_SRC = searchtype,
                   CODE_DEFINITIVE = FALSE,
                   SUGG_SPELLING = NA)
  total <- length(recs)
  pb <- winProgressBar(title = paste0("APHIAID>TAXIZE: via ", searchtype," names"), label=recs[1], min = 0, max = total, width = 300)
    df=df[0,]
    for (i in 1:total) {
      cat(paste0("\t\t",recs[i],"\n"), file = "getTaxaIDs.log", append = TRUE)
      setWinProgressBar(pb, i, title = NULL, label = paste0(recs[i]," (", total-i," left)"))
      this <- tryCatch({
        taxize::get_wormsid(
          query = recs[i],
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
      if (is.null(this)) {
         thisrec= data.frame(joincol = recs[i],
                             CODE = NA,
                             CODE_SVC = 'TAXIZE',
                             CODE_TYPE = 'APHIAID',
                             CODE_SRC = searchtype,
                             CODE_DEFINITIVE = FALSE,
                             SUGG_SPELLING = NA)
      } else {
        tmp=data.frame(this)
        tmp = tmp[,c("ids","multiple_matches","pattern_match")]
        
        thisrec = data.frame(multiple_matches = tmp$multiple_matches,
                          pattern_match = tmp$pattern_match,
                          joincol = recs[i],
                          CODE = tmp$ids,
                          CODE_SRC = searchtype,
                          CODE_SVC = 'TAXIZE',
                          CODE_TYPE = 'APHIAID',
                          CODE_DEFINITIVE = NA,
                          SUGG_SPELLING = NA
        )
        
        thisrec[thisrec$multiple_matches==FALSE & thisrec$pattern_match ==FALSE,"CODE_DEFINITIVE"]<-TRUE
        thisrec[thisrec$multiple_matches==TRUE & thisrec$pattern_match ==FALSE,"CODE_DEFINITIVE"]<-FALSE
        if (is.na(thisrec$CODE_DEFINITIVE))browser()
        thisrec$multiple_matches<-NULL
        thisrec$pattern_match<-NULL
        rm(tmp)
      }
      df = rbind(df,thisrec)
    }
    close(pb)
    df=df[df$joincol!="MMMMMMM",]
    names(df)[names(df) == "joincol"] <- col
    return(df)
  }