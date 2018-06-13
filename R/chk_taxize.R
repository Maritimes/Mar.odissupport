#' @title chk_taxize
#' @description This function sends either the scientific or common name of each
#' species off to taxize
#' @param df - df of species for which the APHIAID is either unknown or not
#' definitive
#' @param field - field in the df containing the value to check against the 
#' service
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom taxize get_wormsid
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_taxize <- function(df = NULL,
           field = NULL,
           searchtype = NULL) {
  total <- nrow(df)
  pb <- winProgressBar(title = paste0("APHIAID>TAXIZE: via ", searchtype," names"), label=df[,field][1], min = 0, max = total, width = 300)
    results=df[0,]
    for (i in 1:total) {
      setWinProgressBar(pb, i, title = NULL, label = #paste( round(i/total*100, 0),"% done")
                          paste0(df[,field][i]," (", total-i," left)"))
      this <- tryCatch({
        taxize::get_wormsid(
          query = df[,field][i],
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
         this= data.frame(ID =df[i,"ID"],
                       CODE = NA,
                       CODE_SVC = 'TAXIZE',
                       CODE_TYPE = 'APHIAID',
                       CODE_SRC = searchtype,
                       CODE_DEFINITIVE = FALSE,
                       SUGG_SPELLING = NA)
      } else {
        tmp=data.frame(this)
        tmp = tmp[,c("ids","multiple_matches","pattern_match")]
        
        this = data.frame(ID = df$ID[i],
                          multiple_matches = tmp$multiple_matches,
                          pattern_match = tmp$pattern_match,
                          CODE = tmp$ids,
                          CODE_SRC = searchtype,
                          CODE_SVC = 'TAXIZE',
                          CODE_TYPE = 'APHIAID',
                          CODE_DEFINITIVE = FALSE,
                          SUGG_SPELLING = NA
        )
        rm(tmp)
        this[this$multiple_matches==FALSE & this$match =="found","CODE_DEFINITIVE"]<-T
        this[this$multiple_matches==TRUE & this$pattern_match ==FALSE,"CODE_DEFINITIVE"]<-T
        this$multiple_matches<-NULL
        this$pattern_match<-NULL
      }
      this = merge(df[,c("ID","SCI_COL_CLN","COMM_COL_CLN")],this, by="ID", all.y=TRUE)
      results = rbind(results,this)
    }
    close(pb)
    return(results)
  }