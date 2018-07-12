#' @title chk_ritis
#' @description This function sends either the scientific or common name of each
#' species off to ritis
#' @param df - df of species for which the TSN is either unknown or not
#' definitive
#' @param field - field in the df containing the value to check against the 
#' service
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom ritis search_scientific
#' @importFrom ritis usage
#' @importFrom ritis search_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_ritis <- function(df = NULL,
                       field = NULL,
                       searchtype = NULL) {
  total <- nrow(df)
  pb <- winProgressBar(title = paste0("TSN>RITIS: via ", searchtype," names"), label=df[,field][1], min = 0, max = total, width = 300)
  results=df[0,]
  
  for (i in 1:total) {
    cat(paste0("\tritis|",searchtype,"|",df[,field][i],"\n"), file = "getTaxaIDs.log", append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(df[,field][i]," (", total-i," left)"))
  if (searchtype == 'scientific') {
      this <- tryCatch({
        as.data.frame(ritis::search_scientific(df[,field][i], wt = "json", raw = FALSE))
      },
      error = function(cond) {
      })
      if(is.null(this)){
        this= data.frame(ID =df[i,"ID"],
                         CODE = NA,
                         CODE_SVC = 'RITIS',
                         CODE_TYPE = 'TSN',
                         CODE_SRC = searchtype,
                         CODE_DEFINITIVE = FALSE,
                         SUGG_SPELLING = NA)
      }else if (nrow(this)<1){
        this= data.frame(ID =df[i,"ID"],
                         CODE = NA,
                         CODE_SVC = 'RITIS',
                         CODE_TYPE = 'TSN',
                         CODE_SRC = searchtype,
                         CODE_DEFINITIVE = FALSE,
                         SUGG_SPELLING = NA)
      }else{
        tmp = data.frame(this[,c("combinedName","tsn")])
        
        this = data.frame(ID =  df$ID[i],
                          CODE = tmp$tsn,
                          CODE_SRC = searchtype,
                          CODE_SVC = 'RITIS',
                          CODE_TYPE = 'TSN',
                          CODE_DEFINITIVE = FALSE,
                          SUGG_SPELLING = tmp$combinedName
        )
        rm(tmp)
        for (j in 1:nrow(this)){
          thisDefCheck <- tryCatch({
            data.frame(ritis::usage(this[j,"CODE"]))
          },
          error = function(cond) {
          })
          if(is.null(thisDefCheck)){
            this[j,"CODE_DEFINITIVE"]<-FALSE
          }else{
            this[j,"CODE_DEFINITIVE"]<-ifelse(thisDefCheck["taxonUsageRating"]=="valid",TRUE,FALSE)
          }
        }
      }
      #if a code is definitive and matches existing spelling, keep it
      if (nrow(this[toupper(this$SUGG_SPELLING) == toupper(df[i,"SCI_COL_CLN"]) & this$CODE_DEFINITIVE == T,])>0){
        this = this[toupper(this$SUGG_SPELLING) == toupper(df[i,"SCI_COL_CLN"]),]
      } else{
        # print("Scientific name doesn't match any valid TSNs perfectly.")
      }
      this=merge(df[,c("ID","SCI_COL_CLN","COMM_COL_CLN")], this, by="ID", all.y = T) 
      results = rbind(results,this)
      rm(this)
  }
  if (searchtype == 'common'){   
      this <- tryCatch({
        as.data.frame(ritis::search_common(df[,field][i], wt = "json", raw = FALSE))
      },
      error = function(cond) {
      })
      if(nrow(this)<1){
        this= data.frame(ID =df[i,"ID"],
                         CODE = NA,
                         CODE_SVC = 'RITIS',
                         CODE_TYPE = 'TSN',
                         CODE_SRC = searchtype,
                         CODE_DEFINITIVE = FALSE,
                         SUGG_SPELLING = NA)
      }else{
        tmp = data.frame(this[,c("commonName","tsn")])
        #sometimes same result returned with different case
        tmp$commonName = toupper(tmp$commonName)
        tmp = unique(tmp)
        this = data.frame(ID =  df$ID[i],
                          CODE = tmp$tsn,
                          CODE_SRC = searchtype,
                          CODE_SVC = 'RITIS',
                          CODE_TYPE = 'TSN',
                          CODE_DEFINITIVE = FALSE,
                          SUGG_SPELLING = tmp$commonName
        )
        rm(tmp)
        for (j in 1:nrow(this)){
          thisDefCheck <- tryCatch({
            data.frame(ritis::usage(this[j,"CODE"]))
          },
          error = function(cond) {
          })
          if(is.null(thisDefCheck)){
            this[j,"CODE_DEFINITIVE"]<-FALSE
          }else{
            this[j,"CODE_DEFINITIVE"]<-ifelse(thisDefCheck["taxonUsageRating"]=="valid",TRUE,FALSE)
          }
        }
      }
      #if a code is definitive and matches existing spelling, keep it
      if (nrow(this[toupper(this$SUGG_SPELLING) == toupper(df[i,"COMM_COL_CLN"]) & this$CODE_DEFINITIVE == T,])>0){
        this = this[toupper(this$SUGG_SPELLING) == toupper(df[i,"COMM_COL_CLN"]),]
      } else{
      }
      this=merge(df[,c("ID","SCI_COL_CLN","COMM_COL_CLN")], this, by="ID", all.y = T)
      results = rbind(results,this)
      rm(this)
  }
  }
  close(pb)
  return(results)
}