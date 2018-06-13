#' @title chk_worrms
#' @description This function sends either the scientific or common name of each
#' species off to worrms
#' @param df - df of species for which the APHIAID is either unknown or not
#' definitive
#' @param field - field in the df containing the value to check against the 
#' service
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom worrms wm_records_name
#' @importFrom worrms wm_records_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_worrms <- function(df = NULL,
           field = NULL,
           searchtype = NULL) {
    total <- nrow(df)
    pb <- winProgressBar(title = paste0("APHIAID>WORRMS: via ", searchtype," names"), label=df[,field][1], min = 0, max = total, width = 300)
    results=df[0,]
    for (i in 1:total) {
      setWinProgressBar(pb, i, title = NULL, label = #paste( round(i/total*100, 0),"% done")
                          paste0(df[,field][i]," (", total-i," left)"))
      if (searchtype == 'scientific') {
        this <- tryCatch({
          worrms::wm_records_name(df[,field][i], fuzzy = F, marine_only = T)
        },
        error = function(cond) {
        })
      } else if (searchtype == 'common'){
        this <- tryCatch({
          worrms::wm_records_common(df[,field][i], fuzzy = F, marine_only = T)
        },
        error = function(cond) {
        })
      }
      if(is.null(this)){
         this= data.frame(ID =df[i,"ID"],
                       CODE = NA,
                       CODE_SVC = 'WORRMS',
                       CODE_TYPE = 'APHIAID',
                       CODE_SRC = searchtype,
                       CODE_DEFINITIVE = FALSE,
                       SUGG_SPELLING = NA)
      } else {
        tmp = data.frame(this[,c("status","valid_AphiaID","valid_name")]) #"unacceptreason"
        
        this = data.frame(ID = df$ID[i],
                          status = tmp$status,
                          CODE = tmp$valid_AphiaID,
                          CODE_SRC = searchtype,
                          CODE_SVC = 'WORRMS',
                          CODE_TYPE = 'APHIAID',
                          CODE_DEFINITIVE = FALSE,
                          SUGG_SPELLING = tmp$valid_name
        )
        rm(tmp)
        this[this$status %in% c("accepted","unaccepted") & !is.na(this$CODE),"CODE_DEFINITIVE"]<-TRUE
        #this$valid_AphiaID<-NULL
        this$status<-NULL
        #this$valid_name<-NULL
      }
      this = merge(df[,c("ID","SCI_COL_CLN","COMM_COL_CLN")],this, by="ID", all.y=TRUE)
      results = rbind(results,this)
    }
    close(pb)
    return(results)
  }
