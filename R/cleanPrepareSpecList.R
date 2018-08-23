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
  
  #there are 2 types of filters- those where:
  # a word or phrase should be removed from a record;
  # the presence of the word or phrase means that the feidl should be set to NA
  #the second type are identified with "*_DROP"
  
  all_Filts <- c(
                #clean these entries out of records when they appear
                #taxonomic identifiers
                "\\bDOMAIN\\b", "\\bKINGDOM\\b","\\bPHYLUM\\b","\\bDIVISION\\b",
                "\\bSUBPHYLUM\\b","\\bSUBDIVISION\\b","\\bCLASS\\b",
                "\\bSUBCLASS\\b","\\bSUPERORDER\\b","\\bORDER\\b",
                "\\bSUBORDER\\b","\\bFAMILY\\b","\\bSUBFAMILY\\b",
                "\\bGENUS\\b","\\bSPECIES\\b",
                #descriptions that don't add to ability to ID
                "(\\(NS\\)|\\bNS)","(-|,|\\s)?UNIDENTIFIED.*",
                "(\\)|^|-|,|\\s)?UNID(EN)?(T?)\\.*","UNID\\. FISH",
                #these remove species abbreviations that I can't do anything with
                "([^']\\b[SP]{1,3}\\.?$)","\\b[a-zA-Z]{1,2}\\.")
  allFilts_Drop <-c(
                    #these just indicate that the field is not of a taxa,
                    #but a modifier describing an aspect of an animal - we don't 
                    #want to try to ID with these
                    "\\bBAIT\\b","\\bDIGESTED\\b","\\bNIDENTIFIED PER\\b",
                    "\\bUNIDENTIFIED SPECIES\\b","\\bUNID (FISH|REMAINS)+\\b", 
                    "\\bREMAINS\\b","\\bSURVEY\\b","\\bFSRS -\\b",
                    "\\bRESERVED\\b","\\bPURSE\\b","\\bINVERTEBRATE\\b",
                    "\\bEGG(S?)\\b","\\bINORGANIC DEBRIS\\b","\\bMIXED\\b",
                    "\\bMUCUS\\b","\\bOPERCULUM\\b","\\bFLUID\\b",
                    "\\bLARVAE\\b","\\s?,?ETC\\.?","\\bAND\\b",
                    "^FISH$","^WATER$","^SHARK$","^SAND$")
  
  comm_Filts_Def <- c(all_Filts, comm_Filts, 
                      #remove the brackets themselves, but retain the  info
                      "\\(", "\\)")
  comm_Filts_Drop <-c(allFilts_Drop)
  
  sci_Filts_Def <- c(all_Filts, sci_Filts, 
                     "\\bOBSOLETE\\b", "\\bUNIDENTIFIED\\b",
                      #anything in brackets
                     "\\(.*?\\)")
  
  sci_Filts_Drop <-c(allFilts_Drop, "\\bWHALE\\b","\\bCRAB\\b", 
                     "\\bLOBSTER\\b", "\\bSHRIMP\\b","\\bIRISH MOSS\\b",
                     "\\bSHARK\\b","\\bCOD WORM\\b","\\bCORALS\\b","\\bSKATE\\b",
                     "\\bFINFISHES\\b","\\bGROUNDFISH\\b","\\bPELAGIC FISH\\b",
                     "\\bSAND TUBE\\b")
  
  if (!is.null(sci_col)) {
    #SCI_COL_CLN populated with sci names
    spec_list$SCI_COL_CLN = toupper(spec_list[, sci_col])
    #drop useless records
    spec_list[grepl(x=spec_list$SCI_COL_CLN,ignore.case = TRUE, pattern = paste(sci_Filts_Drop, collapse = "|")),"SCI_COL_CLN"]<-NA
    #remove the filtered items from the records
    for (i in 1:length(sci_Filts_Def)){
      #remove the strings matched by the filters above
      spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = sci_Filts_Def[i]),"SCI_COL_CLN"]<-
        gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = sci_Filts_Def[i]),"SCI_COL_CLN"],
             pattern = sci_Filts_Def[i],replacement = "") 
      #replace multiple spaces with single space               
      spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = " +"),"SCI_COL_CLN"]<-
        gsub(x = spec_list[grepl(x = spec_list$SCI_COL_CLN,ignore.case = T,pattern = " +"),"SCI_COL_CLN"], pattern = " +",replacement = " ") 
      #remove whitespace
      spec_list[,"SCI_COL_CLN"]<-  trimws(spec_list$SCI_COL_CLN)
    }
    #drop really short records
    spec_list[which(nchar(spec_list$SCI_COL_CLN)<4),"SCI_COL_CLN"]<-NA
  } 
  
  if (!is.null(comm_col)) {
    #COMM_COL_CLN populated with common names
    spec_list$COMM_COL_CLN = toupper(spec_list[, comm_col])
    #drop useless records
    spec_list[grepl(x=spec_list$COMM_COL_CLN,ignore.case = TRUE, pattern = paste(comm_Filts_Drop, collapse = "|")),"COMM_COL_CLN"]<-NA
    #remove the filtered items from the records
    for (i in 1:length(comm_Filts_Def)){
      #remove the strings matched by the filters above
      spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = comm_Filts_Def[i]),"COMM_COL_CLN"]<-
        gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = comm_Filts_Def[i]),"COMM_COL_CLN"],
             pattern = comm_Filts_Def[i],replacement = "") 
      #replace multiple spaces with single space               
      spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = " +"),"COMM_COL_CLN"]<-
        gsub(x = spec_list[grepl(x = spec_list$COMM_COL_CLN,ignore.case = T,pattern = " +"),"COMM_COL_CLN"], pattern = " +",replacement = " ") 
      #remove whitespace
      spec_list[,"COMM_COL_CLN"]<-  trimws(spec_list$COMM_COL_CLN)
    }
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
    #final remove whitespace and drop really short records
    spec_list[which(nchar(spec_list$COMM_COL_CLN)<4),"COMM_COL_CLN"]<-NA
  } 
  #remove Rows that get added if you try to run on non-existent recs
  spec_list =spec_list[!grepl("^NA", rownames(spec_list)),]
  return(spec_list)
}