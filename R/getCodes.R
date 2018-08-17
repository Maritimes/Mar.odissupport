#' @title getCodes
#' @description This function sends names to webservices to try to get codes
#' @param mysteryAPHIAID - dataframe with names for which you want to find AphiaIDs
#' @param definitiveAPHIAID - dataframe with names for which you already have 
#' good AphiaIDs.  The results will be bound to this.
#' @param APHIAIDs_checked - a vector of APHIAIDs that chave already been sent 
#' to the webservices (to avoid duplicate checks)
#' @param mysteryTSN - dataframe with names for which you want to find TSNs
#' @param definitiveTSN - dataframe with names for which you already have 
#' good TSNs.  The results will be bound to this.
#' @param TSNs_checked - a vector of TSNs that chave already been sent 
#' to the webservices (to avoid duplicate checks)
#' @param searchtype - flag indicating whether scientific or common names should 
#' be used checking the services
#' @param thisCode - the particular code that is being looked for (e.g. APHIAID)
#' @param codes - the codes the user has requested (e.g. c("APHIAID,"TSN"))
#' @param masterlist  This is the originally submitted list of species (with an 
#' internally assigned ID)
#' @param logName - this is the name of the logfile in the working directory 
#' that progress should be appended to.
getCodes <- function(mysteryAPHIAID = NULL,
                     definitiveAPHIAID = NULL,
                     APHIAIDs_checked= NULL, 
                     mysteryTSN = NULL,
                     definitiveTSN = NULL,
                     TSNs_checked= NULL,
                     searchtype = NULL,
                     thisCode = NULL, 
                     codes = NULL,
                     masterlist = NULL,
                     logName = logName){
  reviewSpelling <-function(df){
    #this function just drops suggested spelling that match what we already have
    if ("SCI_COL_CLN" %in% names(df) & nrow(df)>0)  df[which(df$SCI_COL_CLN==trimws(gsub(df$SUGG_SPELLING, pattern = "\\(SCIENTIFIC\\)", replacement = ""))),"SUGG_SPELLING"]<-NA
    if ("COMM_COL_CLN" %in% names(df) & nrow(df)>0) df[which(df$COMM_COL_CLN==trimws(gsub(df$SUGG_SPELLING, pattern = "\\(COMMON\\)", replacement = ""))),"SUGG_SPELLING"]<-NA
    return(df)
  }
  
  if (is.null(definitiveAPHIAID)) definitiveAPHIAID=mysteryAPHIAID[0,]
  if (is.null(definitiveTSN)) definitiveTSN=mysteryTSN[0,]
  
  if (thisCode=="APHIAID"){
    dfMystery = mysteryAPHIAID
    dfMysteryOther = mysteryTSN
    dfDefinitive = definitiveAPHIAID
    dfDefinitiveOther = definitiveTSN
  }else{
    dfMystery = mysteryTSN
    dfMysteryOther = mysteryAPHIAID
    dfDefinitive = definitiveTSN
    dfDefinitiveOther = definitiveAPHIAID
  }
  
  start_Timer = Sys.time()
  desc = paste0(thisCode, " (", searchtype, ")")
  cat(paste0(desc, " checks starting\n"), file = logName, append = TRUE)
  
  chkField<-ifelse(searchtype=="scientific","SCI_COL_CLN","COMM_COL_CLN")
  
  #need to adjust in case not  all fields provided
  masterList = masterlist[,names(masterlist) %in% c("ID","SCI_COL_CLN","COMM_COL_CLN")]
  
  if (thisCode=="APHIAID"){
    #TAXIZE
    if (nrow(dfMystery[!is.na(dfMystery[chkField]),])>0){
      forCheckA = dfMystery[!is.na(dfMystery[chkField]),]
      notCheckA = dfMystery[is.na(dfMystery[chkField]),]
      
      if (nrow(forCheckA)>0){
        tmp_sci_taxize = do_taxize(forCheckA,chkField,logName,searchtype)
        defCheck = assignDefinitive(df = tmp_sci_taxize, masterList)
        newdefinitive= defCheck[[1]]
        mysteryAPHIAID= defCheck[[2]]
        rm(tmp_sci_taxize)
      if (nrow(dfDefinitive)<1){
        dfDefinitive<-newdefinitive
      }else{
        dfDefinitive <- unique(rbind(dfDefinitive[dfDefinitive$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
      }
      dfMystery = rbind(mysteryAPHIAID,notCheckA)
    }
    #WORRMS - science
    if (nrow(dfMystery)>0){
      tmp_sci_worrms = do_worrms(dfMystery,chkField,logName,searchtype)
      defCheck = assignDefinitive(df = tmp_sci_worrms, masterList)
      newdefinitive= defCheck[[1]]
      dfMystery= defCheck[[2]]
      rm(tmp_sci_worrms)
      if (nrow(dfDefinitive)<1){
        dfDefinitive<-newdefinitive
      }else{
        dfDefinitive <- unique(rbind(dfDefinitive[dfDefinitive$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
      }
    }
  }else if (thisCode=="TSN"){
    #ritis
    if (nrow(dfMystery[!is.na(dfMystery[chkField]),])>0){
      forCheckT = dfMystery[!is.na(dfMystery[chkField]),]
      notCheckT = dfMystery[is.na(dfMystery[chkField]),]
      if (nrow(forCheckT)>0){
        tmp_sci_ritis = do_ritis(forCheckT,chkField,logName,searchtype)
        defCheck = assignDefinitive(df = tmp_sci_ritis, masterList)
        newdefinitive= defCheck[[1]]
        mysteryTSN= defCheck[[2]]
        rm(tmp_sci_ritis)
        if (nrow(dfDefinitive)<1){
          dfDefinitive<-newdefinitive
        }else{
          dfDefinitive <- unique(rbind(dfDefinitive[dfDefinitive$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        }
      }
      dfMystery = rbind(mysteryTSN, notCheckT)
    }
  }
  #if a mystery code matches the spelling of the sent value, we'll use  that one and drop the others
  if(nrow(dfDefinitive)>0) dfDefinitive<-reviewSpelling(dfDefinitive)
  dfMystery<-reviewSpelling(dfMystery)
  uncertainThis <- dfMystery[is.na(dfMystery$SUGG_SPELLING),] 
  dfMystery <- dfMystery[!(dfMystery$ID %in% uncertainThis$ID),]
  dfMystery<-rbind(uncertainThis,dfMystery)
  rm(uncertainThis)
  
  if("APHIAID" %in% codes & thisCode !="APHIAID"){
    if("TSN" %in% codes){
      recs_W_codeT = dfDefinitive[!is.na(dfDefinitive$CODE) & !(dfDefinitive$ID %in% dfDefinitiveOther$ID),]
      #make sure we don't send the same codes off repeatedly
      recs_W_codeT = recs_W_codeT[!(recs_W_codeT$CODE %in% TSNs_checked),]
      if (nrow(recs_W_codeT)>0){
        start_APHIAIDTSN = Sys.time()
        cat(paste0("\t\tUsing found ",thisCode,"s to search for AphiaID\n"), file = logName, append = TRUE)
        code_APHIAIDRes =   do_worrmsAphiaID(recs_W_codeT$CODE,recs_W_codeT, logName=logName)
        #some field renaming below to prevent a warning
        names(code_APHIAIDRes)[names(code_APHIAIDRes) == 'CODE'] <- 'APHIAID'
        aphiaIDFields = c("ID", "CODE")
        code_APHIAIDRes = merge(recs_W_codeT[,aphiaIDFields],code_APHIAIDRes, by.x = "CODE", by.y="TSN")
        code_APHIAIDRes$CODE<-NULL
        names(code_APHIAIDRes)[names(code_APHIAIDRes) == 'APHIAID'] <- 'CODE'
        defCheck = assignDefinitive(df = code_APHIAIDRes, masterList)
        newdefinitive= defCheck[[1]]
        mysteryOtherFromCode= defCheck[[2]]
        rm(code_APHIAIDRes)
        dfDefinitiveOther <- unique(rbind(dfDefinitiveOther[dfDefinitiveOther$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        #this ensures that the mystery values include the ones that didn't get sent have an aphiaid
        dfMysteryOther = dfMysteryOther[!(dfMysteryOther$ID %in% dfDefinitiveOther$ID),]
        dfMysteryOther = dfMysteryOther[!(dfMysteryOther$ID %in% mysteryOtherFromCode$ID),]
        dfMysteryOther =  rbind(dfMysteryOther, mysteryOtherFromCode)
        
        #if a mystery code matches the spelling of the sent value, we'll use  that one and drop the others
        dfDefinitiveOther<-reviewSpelling(dfDefinitiveOther)
        dfMysteryOther<-reviewSpelling(dfMysteryOther)
        uncertainOther <- dfMysteryOther[is.na(dfMysteryOther$SUGG_SPELLING),]
        dfMysteryOther <- dfMysteryOther[!(dfMysteryOther$ID %in% uncertainOther$ID),]
        dfMysteryOther<-rbind(uncertainOther,dfMysteryOther)
        rm(uncertainOther)
        TSNs_checked = c(TSNs_checked, recs_W_codeT$CODE)
        cat(paste0("\t\tAphiaID (via TSN) checks completed in ",format(.POSIXct(difftime(Sys.time(), start_APHIAIDTSN, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
      }    
    }
  }
  
  if("TSN" %in% codes & thisCode !="TSN"){
    if("APHIAID" %in% codes){
      recs_W_codeA = dfDefinitive[!is.na(dfDefinitive$CODE) & !(dfDefinitive$ID %in% dfDefinitiveOther$ID),]
      #make sure we don't send the same codes off repeatedly
      recs_W_codeA = recs_W_codeA[!(recs_W_codeA$CODE %in% APHIAIDs_checked),]
      if (nrow(recs_W_codeA)>0){
        start_TSNAPHIAID = Sys.time()
        cat(paste0("\t\tUsing found ",thisCode,"s to search for TSNs\n"), file = logName, append = TRUE)
        code_TSNRes = do_worrmsTSN(recs_W_codeA$CODE,recs_W_codeA, logName=logName)
        #some field renaming below to prevent a warning
        names(code_TSNRes)[names(code_TSNRes) == 'CODE'] <- 'TSN'
        tsnFields = c("ID", "CODE")
        code_TSNRes = merge(recs_W_codeA[,tsnFields],code_TSNRes, by.x = "CODE", by.y="APHIAID")
        code_TSNRes$CODE<-NULL
        names(code_TSNRes)[names(code_TSNRes) == 'TSN'] <- 'CODE'
        defCheck = assignDefinitive(df = code_TSNRes, masterList)
        newdefinitive= defCheck[[1]]
        mysteryOtherFromCode= defCheck[[2]]
        rm(code_TSNRes)
        dfDefinitiveOther <- unique(rbind(dfDefinitiveOther[dfDefinitiveOther$CODE_DEFINITIVE %in% TRUE,],newdefinitive))
        #this ensures that the mystery values include the ones that didn't get sent have an aphiaid
        dfMysteryOther = dfMysteryOther[!(dfMysteryOther$ID %in% dfDefinitiveOther$ID),]
        dfMysteryOther = dfMysteryOther[!(dfMysteryOther$ID %in% mysteryOtherFromCode$ID),]
        dfMysteryOther =  rbind(dfMysteryOther, mysteryOtherFromCode)
        
        #if a mystery code matches the spelling of the sent value, we'll use  that one and drop the others
        dfDefinitiveOther<-reviewSpelling(dfDefinitiveOther)
        dfMysteryOther<-reviewSpelling(dfMysteryOther)
        uncertainOther <- dfMysteryOther[is.na(dfMysteryOther$SUGG_SPELLING),] 
        dfMysteryOther <- dfMysteryOther[!(dfMysteryOther$ID %in% uncertainOther$ID),]
        dfMysteryOther<-rbind(uncertainOther,dfMysteryOther)
        rm(uncertainOther) 
        APHIAIDs_checked =  c(APHIAIDs_checked, recs_W_codeA$CODE)
        cat(paste0("\t\tTSN (via ",thisCode,") checks completed in ",format(.POSIXct(difftime(Sys.time(), start_TSNAPHIAID, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
      }
    }
  }
  
  if (thisCode=="APHIAID"){
    mysteryAPHIAID = dfMystery 
    mysteryTSN = dfMysteryOther 
    definitiveAPHIAID = dfDefinitive  
    definitiveTSN = dfDefinitiveOther 
  }else{
    mysteryTSN = dfMystery 
    mysteryAPHIAID = dfMysteryOther 
    definitiveTSN = dfDefinitive 
    definitiveAPHIAID = dfDefinitiveOther 
  }
  
  
  cat(paste0(desc, " checks completed in ",format(.POSIXct(difftime(Sys.time(), start_Timer, units="secs"),tz="GMT"), "%H:%M:%S"),"\n"), file = logName, append = TRUE)
  res = list(mysteryAPHIAID, definitiveAPHIAID, mysteryTSN, definitiveTSN, APHIAIDs_checked, TSNs_checked)
  return(res)
}