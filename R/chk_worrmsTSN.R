#' @title chk_worrmsTSN
#' @description This function sends APHIAIDs to worrms to see if it can find a 
#' corresponding TSN
#' @param df - df of species for which we already have an APHIAID
#' @importFrom worrms wm_external
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_worrmsTSN<-function(df = NULL){
  results=df[0,]
  names(df)[names(df)=="CODE"]= "AphiaID"
  names(df)[names(df)=="CODE_DEFINITIVE"]= "AphiaID_DEF"
  
  #added loop since one record with no results botched the whole call;
  #likely a huge performance hit
  for (i in 1:nrow(df)) {
    this <- tryCatch(
      {
        worrms::wm_external(id = as.integer(df[df$CODE_TYPE == "APHIAID","AphiaID"][i]))
      },
      error=function(cond){
      }
    )
    if(is.null(this)){
      this= data.frame(ID =df[i,"ID"],
                       CODE = NA,
                       CODE_SVC = 'WORRMS',
                       CODE_TYPE = 'TSN',
                       CODE_SRC = 'APHIAID',
                       CODE_DEFINITIVE = FALSE,
                       SUGG_SPELLING = NA)
    }else{
      this = data.frame(CODE = this,
                        ID = df$ID[i],
                        CODE_SVC = 'WORRMS',
                        CODE_TYPE = 'TSN',
                        CODE_SRC = "APHIAID",
                        CODE_DEFINITIVE = FALSE,
                        SUGG_SPELLING = NA
      )
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
    this = merge(df[,c("ID","SCI_COL_CLN","COMM_COL_CLN")],this, by="ID", all.y=TRUE)
    results = rbind(results,this)
  }
  return(results)
}
