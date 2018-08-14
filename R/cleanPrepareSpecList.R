#' @title cleanPrepareSpecList
#' @description This function takes the provided species list and cleans up 
#' common issues with field names
#' @param spec_list - df containing fields you want to check webservices for
#' codes for
#' @param sci_col - name of the field in the df which contains scientific names
#' @param comm_col - name of the field in the df which contains common names
#' @param sci_Filts - regex values that can identify values within the 
#' scientific name column that you want the script to totally ignore
#' @param comm_Filts - regex values that can identify values within the 
#' common name column that you want the script to totally ignore
#' @family speciesCodes
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
cleanPrepareSpecList<-function(spec_list = NULL,
                        sci_col = NULL,
                        comm_col = NULL,
                        sci_Filts = NULL,
                        comm_Filts = NULL
                        ){
  
  spec_list$ID <- seq.int(nrow(spec_list))
  
  if (!is.null(sci_col)) spec_list$SCI_COL_CLN = NA
  if (!is.null(comm_col)) spec_list$COMM_COL_CLN = NA
  
  allFilts <- c("BAIT", "DIGESTED","UNIDENTIFIED PER","UNIDENTIFIED SPECIES",
                "UNID (FISH|REMAINS)+", "REMAINS","SURVEY","FSRS -","RESERVED",
                "PURSE","^FISH( AND| \\,|$)","\\,?\\s?EGG(S?)-?","\\s?LARVAE",
                "INVERTEBRATE","WATER","FLUID","^SAND$","EGGS",
                "INORGANIC DEBRIS","MIXED","MUCUS","OPERCULUM","^SHARK$")
  
  comm_Filts_Def <- c(comm_Filts,"([^']\\b[SP]{1,3}\\.?$)",
                 "([^']\\b[a-zA-Z]{1,2}\\.?$)","SHARK ","^FISH$","/")
  
  sci_Filts_Def <- c(sci_Filts, "WHALE","CETACEAN","CRAB", "LOBSTER","SHRIMP",
                "IRISH MOSS","SHARK","COD WORM","SEA CORALS","SKATE","OBSOLETE",
                "FINFISHES","GROUNDFISH","PELAGIC FISH","\\bAND\\b","SAND TUBE",
                "UNIDENTIFIED")
  
  if (!is.null(sci_col)) {
    #remove whitespace
    spec_list$SCI_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list[, sci_col]))
    #remove recs matching allFilts and sci_filts
    spec_list[grepl(x=spec_list[, sci_col],ignore.case = TRUE, pattern = paste(c(allFilts, sci_Filts_Def), collapse = "|")),"SCI_COL_CLN"]<-NA
    
    #remove bad bits, but retain the rest of the string)
    #SCI names never have brackets - get rid of them, and everything they contain
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\(.*?\\)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\(.*?\\)"),"SCI_COL_CLN"],
           pattern = "\\(.*?\\)",replacement = "") 
    #in case there we nested group, we might have an unmatched, dangling bracket
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\)"),"SCI_COL_CLN"],
           pattern = "\\)",replacement = "") 
    
    #(NS) or NS
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"SCI_COL_CLN"],
           pattern = "(\\(NS\\)|\\bNS)",replacement = "") 
    #SP, SP., SPP and SPP.
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "([^']\\b[SP]{1,3}\\.?$)"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "([^']\\b[SP]{1,3}\\.?$)"),"SCI_COL_CLN"],
           pattern = "([^']\\b[SP]{1,3}\\.?$)",replacement = "") 
    
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\b[a-zA-Z]{1,2}\\."),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = "\\b[a-zA-Z]{1,2}\\."),"SCI_COL_CLN"],
           pattern = "\\b[a-zA-Z]{1,2}\\.",replacement = "") 
    
    #final remove whitespace and drop really short records
    spec_list[,"SCI_COL_CLN"]<-  gsub("(^\\s+)|(\\s+$)", "", spec_list$SCI_COL_CLN)
    spec_list[which(nchar(spec_list$SCI_COL_CLN)<4),"SCI_COL_CLN"]<-NA
    
    #find cases of multiple spaces and replace with single space               
    spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = " +"),"SCI_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = " +"),"SCI_COL_CLN"], pattern = " +",replacement = " ") 
  } 
  if (!is.null(comm_col)) {
    #remove whitespace
    spec_list$COMM_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list[, comm_col]))
    #remove recs matching allFilts and comm_filts
    spec_list[grepl(x=spec_list[, comm_col],ignore.case = TRUE, pattern = paste(c(allFilts, comm_Filts_Def), collapse = "|")),"COMM_COL_CLN"]<-NA
    #remove bad bits, but retain the rest of the string)
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(-|,|\\s)?UNIDENTIFIED.*"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(-|,|\\s)?UNIDENTIFIED.*"),"COMM_COL_CLN"],
           pattern = "(-|,|\\s)?UNIDENTIFIED.*",replacement = "") 
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*"),"COMM_COL_CLN"],
           pattern = "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*",replacement = "")
    #(NS) or NS
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "(\\(NS\\)|\\bNS)"),"COMM_COL_CLN"],
           pattern = "(\\(NS\\)|\\bNS)",replacement = "") 
    
    #stupid "FISH" record
    spec_list[grepl(x = spec_list[,comm_col],ignore.case = T,pattern = "UNID\\. FISH"),"COMM_COL_CLN"]<-NA
    
    #remove "ETC."
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "\\s?,?ETC\\.?"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = "\\s?,?ETC\\.?"),"COMM_COL_CLN"],
           pattern = "\\s?,?ETC\\.?",replacement = "") 
    
    #clunky way of reversing strings so we can search for thing like "Greenland Shark" instead of "Shark, Greenland"
    splitter = spec_list[grepl("(.*),(.*)", spec_list$COMM_COL_CLN),c("ID","COMM_COL_CLN")]
    if (nrow(splitter)>0){
      splitter$one=gsub(",.*$", "", splitter$COMM_COL_CLN)
      splitter$two=gsub("^.*,", "", splitter$COMM_COL_CLN)
      splitter$new=paste0(splitter$two, " ",splitter$one)
      splitter$one<-splitter$two<-splitter$COMM_COL_CLN<-NULL
      spec_list = merge(spec_list,splitter, by = "ID", all.x = T)
      spec_list[!is.na(spec_list$new),"COMM_COL_CLN"]<-spec_list[!is.na(spec_list$new),"new"]
      spec_list$new<-NULL
    }
    
    #find cases of multiple spaces and replace with single space               
    spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = " +"),"COMM_COL_CLN"]<-
      gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = " +"),"COMM_COL_CLN"],
           pattern = " +",replacement = " ") 
    
    #final remove whitespace and drop really short records
    spec_list$COMM_COL_CLN = gsub("(^\\s+)|(\\s+$)", "", toupper(spec_list$COMM_COL_CLN))
    spec_list[which(nchar(spec_list$COMM_COL_CLN)<4),"COMM_COL_CLN"]<-NA
  } 
  #remove Rows that get added if you try to run on non-existent recs
  spec_list =spec_list[!grepl("^NA", rownames(spec_list)),]
  return(spec_list)
}