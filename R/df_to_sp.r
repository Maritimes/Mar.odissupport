#' @title df_to_sp
#' @description This function takes a dataframe and converts in to a spatialpointsdataframe
#' @param df default is \code{NULL}. This is the dataframe you want to spatialize.
#' @param lat.field the default is \code{"LATITUDE"}. the name of the field holding latitude values
#' (in decimal degrees)
#' @param lon.field the default is \code{"LONGITUDE"}.  the name of the field holding longitude
#' values (in decimal degrees)
#' @param the.CRS the default is \code{"+init=epsg:4326"}. This identifies the projection for the input data
#' @return spatialpointsdataframe
#' @family data
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom sp CRS
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
df_to_sp <- function(df = NULL, lat.field="LATITUDE", lon.field="LONGITUDE", the.CRS = '+init=epsg:4326'){
  df.sp = SpatialPointsDataFrame(
    coords = df[, c(lon.field, lat.field)],
    data = df,
    proj4string = CRS(the.CRS)
  )
  return(df.sp)
}
