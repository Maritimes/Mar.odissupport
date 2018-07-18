#' @title chk_ritis
#' @description This function sends either the scientific or common name of each
#' species off to ritis
#' @param recs - vector of taxa names for which the TSN is either unknown or not
#' definitive
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @importFrom ritis search_scientific
#' @importFrom ritis usage
#' @importFrom ritis search_common
#' @importFrom utils winProgressBar
#' @importFrom utils setWinProgressBar
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
chk_ritis <- function(recs = NULL,
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
  pb <- winProgressBar(title = paste0("TSN>RITIS: via ", searchtype," names"), label=recs[1], min = 0, max = total, width = 300)
  results=df[0,]
  for (i in 1:total) {
    cat(paste0("\t\t",recs[i],"\n"), file = "getTaxaIDs.log", append = TRUE)
    setWinProgressBar(pb, i, title = NULL, label = paste0(recs[i]," (", total-i," left)"))
    if (searchtype == 'scientific') {
      this <- tryCatch({
        as.data.frame(ritis::search_scientific(recs[i], wt = "json", raw = FALSE))
      },
      error = function(cond) {
      })
      if(is.null(this)){
        thisrec= data.frame(joincol = trimws(toupper(recs[i])),
                            CODE = NA,
                            CODE_SVC = 'RITIS',
                            CODE_TYPE = 'TSN',
                            CODE_SRC = searchtype,
                            CODE_DEFINITIVE = FALSE,
                            SUGG_SPELLING = NA)
      }else if (nrow(this)<1){
        thisrec= data.frame(joincol = trimws(toupper(recs[i])),
                            CODE = NA,
                            CODE_SVC = 'RITIS',
                            CODE_TYPE = 'TSN',
                            CODE_SRC = searchtype,
                            CODE_DEFINITIVE = FALSE,
                            SUGG_SPELLING = NA)
      }else{
        tmp = data.frame(this[,c("combinedName","tsn")])
        thisrec = data.frame(joincol = trimws(toupper(recs[i])),
                             CODE = tmp$tsn,
                             CODE_SRC = searchtype,
                             CODE_SVC = 'RITIS',
                             CODE_TYPE = 'TSN',
                             CODE_DEFINITIVE = FALSE,
                             SUGG_SPELLING = trimws(toupper(tmp$combinedName))
        )
        #ritis does fuzzy matches which gives too many bad results; this limits 
        #to only those with our spelling
        thisrec = thisrec[thisrec$joincol==thisrec$SUGG_SPELLING,]
        if (nrow(thisrec)>0){
          for (j in 1:nrow(thisrec)){
            thisDefCheck <- tryCatch({
              data.frame(ritis::usage(thisrec[j,"CODE"]))
            },
            error = function(cond) {
            })
            #ritis doesn't tell you how fuzzy the matches are, so the only way to 
            #be certain of the tsn is if 1) the tsn is valid, and 2) that the 
            #returned name is identcal to what was sent.
            if(is.null(thisDefCheck)){
              thisrec[j,"CODE_DEFINITIVE"]<-FALSE
            }else if (thisDefCheck["taxonUsageRating"]=="valid"){
              thisrec[j,"CODE_DEFINITIVE"]<-TRUE
            }
          }
        }else{
          thisrec=df[1,]
        }
        results = rbind(results,thisrec)
        rm(thisrec)
      }
    }
    if (searchtype == 'common'){   
      this <- tryCatch({
        as.data.frame(ritis::search_common(recs[i], wt = "json", raw = FALSE))
      },
      error = function(cond) {
      })
      if(nrow(this)<1){
        thisrec= data.frame(joincol = trimws(toupper(recs[i])),
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
        thisrec = data.frame(joincol = trimws(toupper(recs[i])),
                             CODE = tmp$tsn,
                             CODE_SRC = searchtype,
                             CODE_SVC = 'RITIS',
                             CODE_TYPE = 'TSN',
                             CODE_DEFINITIVE = FALSE,
                             SUGG_SPELLING = trimws(toupper(tmp$commonName))
        )
        #ritis does fuzzy matches which gives too many bad results; this limits 
        #to only those with our spelling
        thisrec = thisrec[thisrec$joincol==thisrec$SUGG_SPELLING,]
        if (nrow(thisrec)>0){
          for (j in 1:nrow(thisrec)){  
            thisDefCheck <- tryCatch({
              data.frame(ritis::usage(thisrec[j,"CODE"]))
            },
            error = function(cond) {
            })
            #ritis doesn't tell you how fuzzy the matches are, so the only way to 
            #be certain of the tsn is if 1) the tsn is valid, and 2) that the 
            #returned name is identcal to what was sent.
            if(is.null(thisDefCheck)){
              thisrec[j,"CODE_DEFINITIVE"]<-FALSE
            }else if (thisDefCheck["taxonUsageRating"]=="valid"){
              thisrec[j,"CODE_DEFINITIVE"]<-TRUE
            }
          }
        }else{
          thisrec=df[1,]
        }
        results = rbind(results,thisrec)
        rm(thisrec)
      }
    }
  }
  close(pb)
  return(results)
}