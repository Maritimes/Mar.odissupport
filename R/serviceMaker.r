#' @title serviceMaker
#' @description This function populates a geodatabase with features extracted from the groundfish database.  For every combination of
#' SEASON/YEAR/SPECIES, a layer is created with all of the null sets
#' @param channel this is an open RODBC connection to Oracle, with access to GROUNDFISH.FGP_TOWS2 and GROUNDFISH.FGP_TOWS_NW2
#' @param write.path This is the path to the geodatabase you will write the data to.
#' @return NULL
#' @family data
#' @importFrom arcgisbinding arc.write
#' @importFrom RODBC sqlQuery
#' @importFrom Mar.utils df_to_sp
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
serviceMaker <- function(channel, write.path){
  #all tows and all catches ever
  tows = sqlQuery(channel, "SELECT * FROM GROUNDFISH.FGP_TOWS2")
  catch = sqlQuery(channel, "SELECT * FROM GROUNDFISH.FGP_TOWS_NW2")
  #find the tows for each discrete survey "series"
  surveys = unique(tows[,c("SURVEYYEAR", "SEASON")])
  for (i in 1:nrow(surveys)){
    #get all the tows for this survey
    surveytows = tows[tows$SURVEYYEAR == surveys$SURVEYYEAR[i] & tows$SEASON == surveys$SEASON[i],]
    #get all catch for this survey
    surveycatch = catch[catch$SURVEYYEAR == surveys$SURVEYYEAR[i] & catch$SEASON == surveys$SEASON[i],]
    if (nrow(surveycatch)<1){
      cat(paste0("\nSkipped ",surveys$SURVEYYEAR[i], surveys$SEASON[i]))
      next
    }
    #get all the species caught on this survey
    surveyspecies = unique(surveycatch[,c("SPECIES")])
    #get rid of 0 values
    surveyspecies = surveyspecies[surveyspecies>0]
    for (s in 1:length(surveyspecies)){
      #if (surveyspecies[s]==0)
      #get the records for catches of this species
      surveycatch_sp = surveycatch[surveycatch$SPECIES == surveyspecies[s],]
      #join the catch records for this species to the tows for this survey
      this_survey = merge(surveytows, surveycatch_sp, all.x = TRUE)
      #cosmetics - populate species-related fields with the values
      this_survey$TAXONOMICNAMEAUTHOR = surveycatch_sp$TAXONOMICNAMEAUTHOR[1]
      this_survey$SCIENTIFICNAME = surveycatch_sp$SCIENTIFICNAME[1]
      this_survey$TAXONOMICSERIALNUMBER = surveycatch_sp$TAXONOMICSERIALNUMBER[1]
      this_survey$SPECIES = surveycatch_sp$SPECIES[1]
      #get rid of NAs -ArcGIS doesn't like them - zeroes instead
      this_survey[is.na(this_survey)] <- 0
      #put the caclulated values for weight and num at top so arcgis sees them
      this_survey<-this_survey[order(-this_survey$TOTALNUMBERSTANDARDIZED, this_survey$TOTALWEIGHTSTANDARDIZED_KG),]
      #write out a csv
      this_survey_season = this_survey$SEASON[1]

      this_survey_sp=Mar.utils::df_to_sp(this_survey, "LATITUDE_DD",  "LONGITUDE_DD" )

      mmm_shapeinfo = list("type"="Point",
                           "hasZ"=FALSE,
                           "hasM" = FALSE,
                           "WKT" = "GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]",
                           "WKID" = 4326)
      thisFile = paste0(this_survey$SEASON[1],"_",this_survey$SURVEYYEAR[1],"_SP",sprintf("%04d",this_survey$SPECIES[1]))
      arc.write(path = paste0(write.path,"/",thisFile), data = this_survey_sp, coords = list(this_survey_sp@data$LONGITUDE_DD, this_survey_sp@data$LATITUDE_DD), shape_info = mmm_shapeinfo)
      cat(paste0("\nDid ", write.path,"/",thisFile))
      }
  }
}
